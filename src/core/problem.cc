// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/mpi.h>
#include <deal.II/base/numbers.h>

#include <prismspf/core/dependencies.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/problem.h>
#include <prismspf/core/simulation_timer.h>
#include <prismspf/core/solution_output.h>
#include <prismspf/core/solve_block.h>
#include <prismspf/core/system_wide.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/triangulation_manager.h>

#include <prismspf/solvers/solver_base.h>
#include <prismspf/solvers/solvers.h>

#include <prismspf/user_inputs/linear_solve_parameters.h>
#include <prismspf/user_inputs/temporal_discretization.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <algorithm>
#include <filesystem>

PRISMS_PF_BEGIN_NAMESPACE

namespace
{
  template <unsigned int dim, unsigned int degree, typename number>
  std::vector<std::shared_ptr<SolverBase<dim, degree, number>>>
  make_solvers(const std::vector<SolveBlock>           &solve_blocks,
               const SolveContext<dim, degree, number> &solve_context)
  {
    // Todo: upgrade to recursive for aux solvers
    std::vector<std::shared_ptr<SolverBase<dim, degree, number>>> solvers;
    solvers.reserve(solve_blocks.size());
    for (const auto &solve_block : solve_blocks)
      {
        switch (solve_block.solve_type)
          {
            case SolveType::Explicit:
              solvers.emplace_back(
                std::make_shared<ExplicitSolver<dim, degree, number>>(solve_block,
                                                                      solve_context));
              break;
            case SolveType::Linear:
              solvers.emplace_back(
                std::make_shared<LinearSolver<dim, degree, number>>(solve_block,
                                                                    solve_context));
              break;
            case SolveType::Newton:
              solvers.emplace_back(
                std::make_shared<NewtonSolver<dim, degree, number>>(solve_block,
                                                                    solve_context));
              break;
            case SolveType::Constant:
              solvers.emplace_back(
                std::make_shared<ConstantSolver<dim, degree, number>>(solve_block,
                                                                      solve_context));
              break;
            default:
              AssertThrow(false, dealii::ExcMessage("Unknown solver type"));
          }
      }
    return solvers;
  }

  std::list<DependencyMap>
  get_all_dependency_sets(const std::vector<SolveBlock> &solve_blocks)
  {
    // Todo: upgrade to recursive for aux solvers
    std::list<DependencyMap> output;
    for (const auto &solve_block : solve_blocks)
      {
        output.push_back(solve_block.dependencies_lhs);
        output.push_back(solve_block.dependencies_rhs);
      }
    return output;
  }

  std::list<SolveBlock>
  get_all_solve_blocks(const std::vector<SolveBlock> &solve_blocks)
  {
    // Todo: upgrade to recursive for aux solvers
    std::list<SolveBlock> output;
    for (const auto &solve_block : solve_blocks)
      {
        output.push_back(solve_block);
      }
    return output;
  }

  /**
   * @brief Check if any solve block uses multigrid.
   * @param solve_blocks The vector of solve blocks to check.
   * @return True if any solve block uses multigrid, false otherwise.
   */
  bool
  has_multigrid(const std::vector<SolveBlock> &solve_blocks)
  {
    return std::any_of(solve_blocks.begin(),
                       solve_blocks.end(),
                       [](const SolveBlock &sb)
                       {
                         return sb.linear_solver_parameters.preconditioner == GMG;
                       });
  }

  template <unsigned int dim, unsigned int degree, typename number>
  std::vector<GroupSolutionHandler<dim, number> *>
  get_solution_managers_from_solvers(
    const std::vector<std::shared_ptr<SolverBase<dim, degree, number>>> &solvers)
  {
    // Todo: upgrade to recursive for aux solvers
    std::vector<GroupSolutionHandler<dim, number> *> solution_managers;
    solution_managers.reserve(solvers.size());
    for (const auto &solver : solvers)
      {
        solution_managers.push_back(&(solver->get_solution_manager()));
      }
    return solution_managers;
  }
} // namespace

// *1 Big TODO: Make these classes default-constructible, then use their `init()`
// functions. Also, it may be wise to make SolveContext the owner of all these as well
// rather than Problem. I just want to get this working, so i'm working around this for
// now. TriangulationManager, DoFManager, ConstraintManager, solvers, SolutionIndexer

template <unsigned int dim, unsigned int degree, typename number>
Problem<dim, degree, number>::Problem(
  const std::vector<FieldAttributes>   &_field_attributes,
  const std::vector<SolveBlock>        &_solve_blocks,
  const UserInputParameters<dim>       &_user_inputs,
  PhaseFieldTools<dim>                 &_pf_tools,
  PDEOperatorBase<dim, degree, number> &_pde_operator)
  : field_attributes(_field_attributes)
  , solve_blocks(validate_solve_blocks(_solve_blocks, _field_attributes))
  , user_inputs_ptr(&_user_inputs)
  , pf_tools(&_pf_tools)
  , solve_context(field_attributes,
                  _user_inputs,
                  triangulation_manager,
                  dof_manager,
                  constraint_manager,
                  solution_indexer,
                  _pde_operator)
  , grid_refiner(solve_context)
{
  std::map<unsigned int, LinearSolverParameters> linear_solver_parameters_copy =
    _user_inputs.linear_solve_parameters.linear_solvers;
  std::map<unsigned int, NonlinearSolverParameters> nonlinear_solver_parameters_copy =
    _user_inputs.nonlinear_solve_parameters.newton_solvers;
  for (auto &solve_block : solve_blocks)
    {
      solve_block.validate();
      if (const auto &lin_param_it = linear_solver_parameters_copy.find(solve_block.id);
          lin_param_it != linear_solver_parameters_copy.end())
        {
          ConditionalOStreams::pout_base()
            << "Overriding linear solver parameters for solve block " << solve_block.id
            << " with user input parameters.\n";
          solve_block.linear_solver_parameters = lin_param_it->second;
          linear_solver_parameters_copy.erase(lin_param_it);
        }
      if (const auto &nonlin_param_it =
            nonlinear_solver_parameters_copy.find(solve_block.id);
          nonlin_param_it != nonlinear_solver_parameters_copy.end())
        {
          ConditionalOStreams::pout_base()
            << "Overriding newton solver parameters for solve block " << solve_block.id
            << " with user input parameters.\n";
          solve_block.nonlinear_solver_parameters = nonlin_param_it->second;
          nonlinear_solver_parameters_copy.erase(nonlin_param_it);
        }
    }
  for (const auto &remaining_lin_params : linear_solver_parameters_copy)
    {
      ConditionalOStreams::pout_base()
        << "Warning: Linear solver parameters provided by user inputs for solve "
           "block "
        << remaining_lin_params.first
        << " which does not exist in the solve blocks. These parameters will be "
           "ignored.\n";
    }
  for (const auto &remaining_nonlin_params : nonlinear_solver_parameters_copy)
    {
      ConditionalOStreams::pout_base()
        << "Warning: Nonlinear solver parameters provided by user inputs for solve "
           "block "
        << remaining_nonlin_params.first
        << " which does not exist in the solve blocks. These parameters will be "
           "ignored.\n";
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
Problem<dim, degree, number>::init_system()
{
  const UserInputParameters<dim> &user_inputs = *user_inputs_ptr;
  const unsigned int n_proc = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  const unsigned int n_vect_doubles = dealii::VectorizedArray<number>::size();
  const unsigned int n_vect_bits    = 8 * sizeof(number) * n_vect_doubles;

  ConditionalOStreams::pout_base() << "number of processes: " << n_proc << "\n"
                                   << std::flush;
  ConditionalOStreams::pout_base()
    << "vectorization over " << n_vect_doubles << " doubles = " << n_vect_bits
    << " bits (" << dealii::Utilities::System::get_current_vectorization_level() << ')'
    << "\n"
    << std::flush;

  bool use_mg = has_multigrid(solve_blocks);

  // Create the mesh
  ConditionalOStreams::pout_base() << "Creating triangulation...\n" << std::flush;
  Timer::start_section("Generate mesh");
  triangulation_manager.generate_mesh(user_inputs.spatial_discretization);
  if (use_mg)
    {
      triangulation_manager.init_mg();
    }
  Timer::end_section("Generate mesh");

  // Create the dof handlers.
  ConditionalOStreams::pout_base() << "Creating DoFHandlers...\n" << std::flush;
  Timer::start_section("reinitialize DoFHandlers");
  dof_manager.reinit(triangulation_manager, use_mg);
  dof_manager.reinit_mapping(field_attributes);
  Timer::end_section("reinitialize DoFHandlers");

  // Create the constraints
  // See *1
  ConditionalOStreams::pout_base() << "Creating constraints...\n" << std::flush;
  Timer::start_section("Create constraints");
  constraint_manager.init(user_inputs.boundary_parameters,
                          user_inputs.spatial_discretization,
                          dof_manager,
                          solve_context.get_pde_operator());
  constraint_manager.reinit(field_attributes);
  Timer::end_section("Create constraints");

  Timer::start_section("Initialize MatrixFree");
  solve_context.get_matrix_free_manager().reinit(solve_context.get_dof_manager(),
                                                 solve_context.get_constraint_manager());
  Timer::end_section("Initialize MatrixFree");

  // InvM
  Timer::start_section("Initialize InvM");
  solve_context.get_invm_manager().reinit(solve_context.get_matrix_free_manager());
  solve_context.get_invm_manager().compute_invm();
  Timer::end_section("Initialize InvM");

  // Initialize the solvers
  Timer::start_section("Initialize Solvers");
  solvers = make_solvers(solve_blocks, solve_context);
  solution_indexer.init(field_attributes.size(),
                        get_solution_managers_from_solvers(solvers));
  for (auto &solver : solvers)
    {
      solver->init(get_all_solve_blocks(solve_blocks));
    }
  Timer::end_section("Initialize Solvers");

  // Update the ghosts
  Timer::start_section("Update ghosts");
  for (auto &solver : solvers)
    {
      solver->update_ghosts();
    }
  Timer::end_section("Update ghosts");

  // Perform the initial grid refinement. For this one, we have to do a loop to sufficient
  // coarsen cells to the minimum level
  ConditionalOStreams::pout_base() << "initializing grid refiner..." << std::flush;
  grid_refiner.add_refinement_marker(
    std::make_shared<NucleusRefinementFunction<dim>>(user_inputs.nucleation_parameters,
                                                     pf_tools->nuclei_list));
}

template <unsigned int dim, unsigned int degree, typename number>
void
Problem<dim, degree, number>::solve()
{
  Timer::start_section("Problem Solve");
  // Print a warning if running in DEBUG mode
  ConditionalOStreams::pout_verbose()
    << "\n\n"
    << "================================================\n"
       "  Warning: running in DEBUG mode \n"
    << "================================================\n\n\n"
    << std::flush;
  ConditionalOStreams::pout_summary()
    << "================================================\n"
       "  Initialization\n"
    << "================================================\n"
    << std::flush;

  Timer::start_section("Initialization");
  init_system();
  Timer::end_section("Initialization");

  ConditionalOStreams::pout_base() << "\nSolving...\n\n" << std::flush;

  ConditionalOStreams::pout_summary()
    << "================================================\n"
       "  Solve\n"
    << "================================================\n"
    << std::flush;

  const UserInputParameters<dim> &user_inputs = *user_inputs_ptr;
  const TemporalDiscretization   &time_info   = user_inputs.temporal_discretization;
  SimulationTimer                &sim_timer   = solve_context.get_simulation_timer();
  // Main time-stepping loop
  int exit_status = 0;
  while (sim_timer.get_increment() <= time_info.num_increments && exit_status == 0)
    {
      // Solve a single increment
      // Includes nucleation, refinement, constraints, solve, output, and update
      exit_status = solve_increment(sim_timer);
      // Update time
      sim_timer.increment();
    }

  // Print summary of nuclei seeded during the simulation
  ConditionalOStreams::pout_summary()
    << "================================================\n"
       "  Nuclei Seeded\n"
    << "================================================\n"
    << std::to_string(pf_tools->nuclei_list.size()) << " total nuclei seeded.\n";
  for (const Nucleus<dim> nucleus : pf_tools->nuclei_list)
    {
      ConditionalOStreams::pout_summary() << nucleus << "\n";
    }
  ConditionalOStreams::pout_summary() << "\n" << std::flush;

  Timer::end_section("Problem Solve");
  // Print timer summary
  Timer::print_summary();

  // Throw exception if we exitied for a bad reason
  switch (exit_status)
    {
      case 0: // normal
      case 1: // exit early as normal behavior
        break;
      case 2: // exit early because NaN
        if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
          {
            throw ExcNaN("Exiting early.\n");
          }
        break;
      case 3: // exit triggered by user
        ConditionalOStreams::pout_base() << "\nExiting triggered by user.\n";
        break;
      default:
        break;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
int
Problem<dim, degree, number>::solve_increment(SimulationTimer &sim_timer)
{
  int                             exit_status  = 0;
  bool                            force_output = false;
  const UserInputParameters<dim> &user_inputs  = *user_inputs_ptr;
  unsigned int                    increment    = sim_timer.get_increment();
  bool is_output_increment = user_inputs.output_parameters.should_output(increment);
  bool is_nucleation_increment =
    user_inputs.nucleation_parameters.should_attempt_nucleation(increment);

  // Update the time-dependent constraints
  Timer::start_section("Update time-dependent constraints");
  // TODO: Loop over levels, pass in current time
  constraint_manager.update_time_dependent_constraints(field_attributes);
  Timer::end_section("Update time-dependent constraints");

  // Solve a single increment
  Timer::start_section("Solvers");
  for (auto &solver : solvers)
    {
      SolveTiming solve_timing = solver->get_solve_block().solve_timing;
      if ((solve_timing == PostProcess && !is_output_increment) ||
          (solve_timing == NucleationRate &&
           !(is_nucleation_increment || is_output_increment)))
        {
          continue;
        }
      solve_context.get_pde_operator().pre_solve_block(solve_context,
                                                       solver->get_solve_block().id);
      solver->solve();
      solver->update_ghosts();
      solve_context.get_pde_operator().post_solve_block(solve_context,
                                                        solver->get_solve_block().id);
    }
  Timer::end_section("Solvers");

  // Check for NaN. This isn't an exhaustive search. Just a quick check on specific
  // values.
  for (unsigned int field_index = 0;
       field_index < solve_context.get_field_attributes().size();
       ++field_index)
    {
      if (dealii::Utilities::MPI::logical_or(!dealii::numbers::is_finite(
                                               solve_context.get_solution_indexer()
                                                 .get_solution_vector(field_index)
                                                 .local_element(0)),
                                             MPI_COMM_WORLD))
        {
          exit_status  = 2;
          force_output = true;
          break;
        }
    }

  // Check for user triggered stop
  if (dealii::Utilities::MPI::logical_or(solve_context.get_pde_operator().get_user_stop(),
                                         MPI_COMM_WORLD))
    {
      exit_status  = 3;
      force_output = true;
    }

  // Check for stochastic nucleation.
  bool any_nucleation_occurred = false;
  if (is_nucleation_increment)
    {
      Timer::start_section("Check for nucleation");
      any_nucleation_occurred =
        NucleationManager<dim, degree, number>::attempt_nucleation(solve_context,
                                                                   pf_tools->nuclei_list);
      Timer::end_section("Check for nucleation");
    }

  // Perform grid refinement if necessary
  if (user_inputs.spatial_discretization.has_adaptivity && increment == 0)
    {
      Timer::start_section("Grid refinement");
      grid_refiner.do_initial_refinement(solvers);
      Timer::end_section("Grid refinement");
    }
  else if (user_inputs.spatial_discretization.has_adaptivity &&
           (user_inputs.spatial_discretization.should_refine_mesh(increment) ||
            any_nucleation_occurred))
    {
      // Perform grid refinement
      ConditionalOStreams::pout_base()
        << "[Increment " << sim_timer.get_increment() << "] : Grid Refinement\n";
      Timer::start_section("Grid refinement");
      grid_refiner.do_adaptive_refinement(solvers);
      Timer::end_section("Grid refinement");
      ConditionalOStreams::pout_base() << "\n" << std::flush;
    }

  // Output results if needed
  if (is_output_increment || force_output)
    {
      std::filesystem::path output_prefix =
        std::filesystem::path(user_inputs.output_parameters.directory) /
        user_inputs.output_parameters.file_name;
      std::filesystem::path output_path = output_prefix;
      output_path.remove_filename();
      std::filesystem::create_directories(output_path);
      Timer::start_section("Output");
      SolutionOutput<dim, degree, number>(field_attributes,
                                          solve_context.get_solution_indexer(),
                                          sim_timer,
                                          dof_manager,
                                          output_prefix,
                                          user_inputs);

      // Print the l2-norms and integrals of each solution
      ConditionalOStreams::pout_base()
        << "Iteration: " << sim_timer.get_increment() << "\n";
      for (unsigned int index = 0; index < field_attributes.size(); ++index)
        {
          const auto &solution =
            solve_context.get_solution_indexer().get_solution_vector(index);
          ConditionalOStreams::pout_base()
            << " Solution index " << index << " l2-norm: " << solution.l2_norm()
            << " integrated value: ";

          const auto tensor_rank = field_attributes[index].field_type;

          if (tensor_rank == TensorRank::Vector)
            {
              ConditionalOStreams::pout_base()
                << Integrator<dim, degree, number>::template integrate<1>(
                     dof_manager.get_field_dof_handler(index),
                     solution)
                << "\n";
            }
          else if (tensor_rank == TensorRank::Scalar)
            {
              // This is equivalent to integration
              ConditionalOStreams::pout_base()
                << solution * solve_context.get_invm_manager().get_jxw(TensorRank::Scalar)
                << "\n";
            }
        }
      ConditionalOStreams::pout_base() << "\n" << std::flush;
      Timer::end_section("Output");
    }

  // Update the field labels in preparation for next increment (c_n -> c_n-1)
  for (auto &solver : solvers)
    {
      solver->update();
    }
  return exit_status;
}

template <unsigned int dim, unsigned int degree, typename number>
const SolveContext<dim, degree, number> &
Problem<dim, degree, number>::get_solve_context() const
{
  return solve_context;
}

#include "core/problem.inst"

PRISMS_PF_END_NAMESPACE
