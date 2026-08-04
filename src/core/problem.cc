// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
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
#include <prismspf/core/triangulation_manager.h>

#include <prismspf/solvers/solver_base.h>
#include <prismspf/solvers/solvers.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/timer.h>

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
                       [](const SolveBlock &solve_block)
                       {
                         return solve_block.linear_solver_parameters.preconditioner ==
                                GMG;
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
  // Override boundary condition parameters if they are specified in user inputs
  std::unordered_map<std::string, BoundaryConditionSet> boundary_condition_list =
    _user_inputs.boundary_parameters.boundary_condition_list;
  for (auto &field : field_attributes)
    {
      if (const auto &bc_it = boundary_condition_list.find(field.name);
          bc_it != boundary_condition_list.end())
        {
          Logger::instance() << LogFormatter::info(
                                  "Overriding boundary condition parameters for field " +
                                  field.name + " with user input parameters")
                             << std::endl;
          field.boundary_conditions = bc_it->second;
          boundary_condition_list.erase(bc_it);
        }
    }
  for (const auto &remaining_bc_params : boundary_condition_list)
    {
      Logger::instance()
        << LogFormatter::warning(
             "Boundary condition parameters provided by user inputs for field " +
             remaining_bc_params.first + " will be ignored!")
        << std::endl;
    }

  // Override solver parameters if they are specified in user inputs
  std::map<unsigned int, LinearSolverParameters> linear_solver_parameters_copy =
    _user_inputs.linear_solve_parameters.linear_solvers;
  std::map<unsigned int, NonlinearSolverParameters> nonlinear_solver_parameters_copy =
    _user_inputs.nonlinear_solve_parameters.nonlinear_solvers;
  for (auto &solve_block : solve_blocks)
    {
      solve_block.validate();
      if (const auto &lin_param_it = linear_solver_parameters_copy.find(solve_block.id);
          lin_param_it != linear_solver_parameters_copy.end())
        {
          Logger::instance() << LogFormatter::info(
                                  "Overriding linear solver parameters for solve block " +
                                  std::to_string(solve_block.id) +
                                  " with user input parameters")
                             << std::endl;
          solve_block.linear_solver_parameters = lin_param_it->second;
          linear_solver_parameters_copy.erase(lin_param_it);
        }
      if (const auto &nonlin_param_it =
            nonlinear_solver_parameters_copy.find(solve_block.id);
          nonlin_param_it != nonlinear_solver_parameters_copy.end())
        {
          Logger::instance()
            << LogFormatter::info(
                 "Overriding nonlinear solver parameters for solve block " +
                 std::to_string(solve_block.id) + " with user input parameters")
            << std::endl;
          solve_block.nonlinear_solver_parameters = nonlin_param_it->second;
          nonlinear_solver_parameters_copy.erase(nonlin_param_it);
        }
    }
  for (const auto &remaining_lin_params : linear_solver_parameters_copy)
    {
      Logger::instance()
        << LogFormatter::warning(
             "Linear solver parameters provided by user inputs for solve block " +
             std::to_string(remaining_lin_params.first) + " will be ignored!")
        << std::endl;
    }
  for (const auto &remaining_nonlin_params : nonlinear_solver_parameters_copy)
    {
      Logger::instance()
        << LogFormatter::warning(
             "Nonlinear solver parameters provided by user inputs for solve block " +
             std::to_string(remaining_nonlin_params.first) + " will be ignored!")
        << std::endl;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
Problem<dim, degree, number>::init_system()
{
  Timer::Scope problem("Init");

  // Print some basic initialization stuff
  {
    const unsigned int n_proc = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    const unsigned int n_vect_doubles = dealii::VectorizedArray<number>::size();
    const unsigned int n_vect_bits    = 8 * sizeof(number) * n_vect_doubles;

    Logger::instance() << LogFormatter::section("Initialization") << std::endl;
#ifdef DEBUG
    Logger::instance() << LogFormatter::warning("Running in debug!") << std::endl;
#endif
    Logger::instance() << LogFormatter::info("Number of processes " +
                                             std::to_string(n_proc))
                       << std::endl;
    Logger::instance() << LogFormatter::info(
                            "Vectorization over " + std::to_string(n_vect_doubles) +
                            " doubles = " + std::to_string(n_vect_bits) + " bits (" +
                            dealii::Utilities::System::get_current_vectorization_level() +
                            ")")
                       << std::endl
                       << std::endl;
  }

  const UserInputParameters<dim> &user_inputs = *user_inputs_ptr;
  bool                            use_mg      = has_multigrid(solve_blocks);

  // Create the mesh
  Logger::instance() << "Generating mesh...\n" << std::flush;
  {
    // TODO: Consolidate these two methods
    triangulation_manager.generate_mesh(user_inputs.spatial_discretization);
    if (use_mg)
      {
        triangulation_manager.init_mg();
      }
  }
  Logger::instance() << LogFormatter::success("Mesh generation succeeded") << std::endl;

  // Create the dof handlers.
  Logger::instance() << "Creating DoFHandlers...\n" << std::flush;
  {
    dof_manager.reinit(triangulation_manager, use_mg);
    dof_manager.reinit_mapping(field_attributes);
  }
  Logger::instance() << LogFormatter::success("DoFHandler creation succeeded")
                     << std::endl;

  // Create the constraints
  // See *1
  Logger::instance() << "Creating constraints...\n" << std::flush;
  {
    constraint_manager.init(user_inputs.boundary_parameters,
                            user_inputs.spatial_discretization,
                            dof_manager,
                            solve_context.get_pde_operator(),
                            solve_context.get_simulation_timer());
    constraint_manager.reinit(field_attributes);
  }
  Logger::instance() << LogFormatter::success("Constraint creation succeeded")
                     << std::endl;

  // Create MatrixFree
  Logger::instance() << "Creating MatrixFree...\n" << std::flush;
  {
    solve_context.get_matrix_free_manager()
      .reinit(solve_context.get_dof_manager(), solve_context.get_constraint_manager());
  }
  Logger::instance() << LogFormatter::success("MatrixFree creation succeeded")
                     << std::endl;

  // Create InvM
  Logger::instance() << "Creating InvM...\n" << std::flush;
  {
    solve_context.get_invm_manager().reinit(solve_context.get_matrix_free_manager());
    solve_context.get_invm_manager().compute_invm();
  }
  Logger::instance() << LogFormatter::success("InvM creation succeeded") << std::endl;

  // Initialize the solvers
  Logger::instance() << "Initializing solvers...\n" << std::flush;
  {
    solvers = make_solvers(solve_blocks, solve_context);
    solution_indexer.init(field_attributes.size(),
                          get_solution_managers_from_solvers(solvers));
    for (auto &solver : solvers)
      {
        solver->init(get_all_solve_blocks(solve_blocks));
      }

    // Update the ghosts
    // TODO: Log this?
    {
      for (auto &solver : solvers)
        {
          solver->update_ghosts();
        }
    }
  }
  Logger::instance() << LogFormatter::success("Solve initialization succeeded")
                     << std::endl;

  // Attach refinement function for nucleation
  // TODO: Log this?
  grid_refiner.add_refinement_marker(
    std::make_shared<NucleusRefinementFunction<dim>>(user_inputs.nucleation_parameters,
                                                     pf_tools->nuclei_list));
}

template <unsigned int dim, unsigned int degree, typename number>
void
Problem<dim, degree, number>::solve()
{
  int exit_status = 0;
  {
    Timer::Scope problem("Problem");

    const UserInputParameters<dim> &user_inputs = *user_inputs_ptr;
    const TemporalDiscretization   &time_info   = user_inputs.temporal_discretization;
    SimulationTimer                &sim_timer   = solve_context.get_simulation_timer();

    // TODO: Remove these asserts as the features/bugs are fixed
    AssertThrow(!user_inputs.spatial_discretization.has_adaptivity || dim != 1,
                dealii::ExcMessage(
                  "AMR cannot be enable for 1D in deal.II 9.7.0 and below."));
    AssertThrow(
      !user_inputs.spatial_discretization.has_adaptivity || !has_multigrid(solve_blocks),
      dealii::ExcMessage(
        "AMR cannot be enabled when using multigrid preconditioners currently."));

    init_system();

    // Main time-stepping loop
    Logger::instance() << LogFormatter::section("Time-stepping Loop") << std::endl;
    while (sim_timer.get_increment() <= time_info.n_increments && exit_status == 0)
      {
        // Solve a single increment
        // Includes nucleation, refinement, constraints, solve, output, and update
        exit_status = solve_increment(sim_timer);
        // Update time
        sim_timer.increment();
      }

    Logger::instance() << LogFormatter::section("Post-simulation Details") << std::endl;
    {
      // Print summary of nuclei seeded during the simulation
      Logger::instance() << LogFormatter::subsection("Nuclei Seeded") << std::endl;
      {
        Logger::instance() << "Total nuclei seeded " << pf_tools->nuclei_list.size()
                           << std::endl;
        // TODO: This only belongs in the log file
        Logger::IndentScope nuclei_scope;
        for (const Nucleus<dim> &nucleus : pf_tools->nuclei_list)
          {
            Logger::instance() << nucleus << std::endl;
          }
      }
    }
  }
  Timer::print_summary();

  // Throw exception if we exitied for a bad reason
  switch (exit_status)
    {
      case 0: // normal
      case 1: // exit early as normal behavior
        break;
      case 2: // exit early because NaN
              // TODO: Use assertion to throw
        if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
          {
            throw ExcNaN("Exiting early.\n");
          }
        break;
      case 3: // exit triggered by user
        Logger::reset_indent();
        Logger::instance() << LogFormatter::error("Exit triggered by user!") << std::endl;
        break;
      default:
        break;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
int
Problem<dim, degree, number>::solve_increment(SimulationTimer &sim_timer)
{
  Timer::Scope problem("Solve Increment");

  int                             exit_status  = 0;
  bool                            force_output = false;
  const UserInputParameters<dim> &user_inputs  = *user_inputs_ptr;
  unsigned int                    increment    = sim_timer.get_increment();
  bool is_output_increment = user_inputs.output_parameters.should_output(increment);
  bool is_nucleation_increment =
    user_inputs.nucleation_parameters.should_attempt_nucleation(increment);

  // Update the time-dependent constraints
  // TODO: Loop over levels, pass in current time
  constraint_manager.update_time_dependent_constraints(field_attributes);

  // Solve a single increment
  {
    Timer::Scope solve_scope("Solve Blocks");
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
  }

  // Check for NaN. This isn't an exhaustive search. Just a quick check on specific
  // values.
  for (unsigned int field_index = 0;
       field_index < solve_context.get_field_attributes().size();
       ++field_index)
    {
      const auto &vec =
        solve_context.get_solution_indexer().get_solution_vector(field_index);
      bool not_finite =
        vec.locally_owned_size() > 0 && !dealii::numbers::is_finite(vec.local_element(0));
      if (dealii::Utilities::MPI::logical_or(not_finite, MPI_COMM_WORLD))
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
      any_nucleation_occurred =
        NucleationManager<dim, degree, number>::attempt_nucleation(solve_context,
                                                                   pf_tools->nuclei_list);
    }

  // Perform grid refinement if necessary
  if (user_inputs.spatial_discretization.has_adaptivity && increment == 0)
    {
      grid_refiner.do_initial_refinement(solvers);
    }
  else if (user_inputs.spatial_discretization.has_adaptivity &&
           (user_inputs.spatial_discretization.should_refine_mesh(increment) ||
            any_nucleation_occurred))
    {
      // Perform grid refinement
      grid_refiner.do_adaptive_refinement(solvers);
    }

  // Output results if needed
  if (is_output_increment || force_output)
    {
      std::filesystem::path output_prefix =
        std::filesystem::path(user_inputs.output_parameters.folder) /
        user_inputs.output_parameters.file_name;
      std::filesystem::path output_path = output_prefix;
      output_path.remove_filename();
      std::filesystem::create_directories(output_path);
      SolutionOutput<dim, degree, number>(field_attributes,
                                          solve_context.get_solution_indexer(),
                                          sim_timer,
                                          dof_manager,
                                          output_prefix,
                                          user_inputs);

      // Print the l2-norms and integrals of each solution
      Logger::instance() << "Iteration " << sim_timer.get_increment() << " Time "
                         << sim_timer.get_time() << std::endl;
      for (unsigned int index = 0; index < field_attributes.size(); ++index)
        {
          Logger::IndentScope output_scope;

          const auto &solution =
            solve_context.get_solution_indexer().get_solution_vector(index);

          Logger::instance() << "Field index " << index << " name "
                             << field_attributes[index].name << " l2-norm "
                             << solution.l2_norm() << " integrated value ";

          const auto tensor_rank = field_attributes[index].field_type;

          if (tensor_rank == TensorRank::Vector)
            {
              Logger::instance()
                << Integrator<dim, degree, number>::template integrate<1>(
                     dof_manager.get_field_dof_handler(index),
                     solution)
                << std::endl;
            }
          else if (tensor_rank == TensorRank::Scalar)
            {
              // This is equivalent to integration
              Logger::instance()
                << solution * solve_context.get_invm_manager().get_jxw(TensorRank::Scalar)
                << std::endl;
            }
        }
      Logger::instance() << std::endl;
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
