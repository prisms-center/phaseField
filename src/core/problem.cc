// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/problem.h>
#include <prismspf/core/simulation_timer.h>
#include <prismspf/core/solution_output.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/system_wide.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/triangulation_manager.h>

#include <prismspf/solvers/solvers.h>

#include <prismspf/user_inputs/temporal_discretization.h>
#include <prismspf/user_inputs/user_input_parameters.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
std::vector<GroupSolutionHandler<dim, number> *>
get_solution_managers_from_solvers(
  const std::vector<std::shared_ptr<GroupSolverBase<dim, degree, number>>> &solvers)
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

template <unsigned int dim, unsigned int degree, typename number>
Problem<dim, degree, number>::Problem(
  const std::vector<FieldAttributes>                          &_field_attributes,
  const std::vector<SolveGroup>                               &_solve_groups,
  const UserInputParameters<dim>                              &_user_inputs,
  PhaseFieldTools<dim>                                        &_pf_tools,
  const std::shared_ptr<PDEOperatorBase<dim, degree, number>> &_pde_operator)
  : field_attributes(_field_attributes)
  , solve_groups(_solve_groups)
  , user_inputs_ptr(&_user_inputs)
  , pf_tools(&_pf_tools)
  , triangulation_manager(_user_inputs.get_spatial_discretization(), false)
  , dof_manager(field_attributes, triangulation_manager)
  , constraint_manager(field_attributes, dof_manager, _pde_operator.get())
  , solvers(solve_groups.size(), nullptr)
  , solution_indexer(field_attributes.size(), get_solution_managers_from_solvers(solvers))
  , solve_context(field_attributes,
                  _user_inputs,
                  triangulation_manager,
                  dof_manager,
                  constraint_manager,
                  solution_indexer,
                  _pde_operator)
  , grid_refiner(solve_context)
{}

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

  // Create the mesh
  ConditionalOStreams::pout_base() << "creating triangulation...\n" << std::flush;
  Timer::start_section("Generate mesh");
  triangulation_manager.generate_mesh(user_inputs.get_spatial_discretization());
  Timer::end_section("Generate mesh");

  // Create the dof handlers.
  ConditionalOStreams::pout_base() << "creating DoFHandlers...\n" << std::flush;
  Timer::start_section("reinitialize DoFHandlers");
  dof_manager.init(triangulation_manager);
  Timer::end_section("reinitialize DoFHandlers");

  // Create the constraints
  ConditionalOStreams::pout_base() << "creating constraints...\n" << std::flush;
  Timer::start_section("Create constraints");
  for (const auto &solve_group : solve_groups)
    {
      (void) solve_group;
      // TODO: Loop over levels, pass in current time
      // constraint_manager.make_constraints(dof_manager.get_field_dof_handlers(
      //                                      solve_group.field_indices));
    }
  Timer::end_section("Create constraints");

  // Initialize the solvers
  Timer::start_section("Initialize Solvers");
  solvers.reserve(solve_groups.size());
  for (const auto &solve_group : solve_groups)
    {
      std::shared_ptr<GroupSolverBase<dim, degree, number>> solver;
      switch (solve_group.pde_type)
        {
          case PDEType::Explicit:
            solver = std::make_shared<ExplicitSolver<dim, degree, number>>(solve_group,
                                                                           solve_context);
            break;
          case PDEType::Linear:
            solver = std::make_shared<LinearSolver<dim, degree, number>>(solve_group,
                                                                         solve_context);
            break;
          case PDEType::Newton:
            solver = std::make_shared<NewtonSolver<dim, degree, number>>(solve_group,
                                                                         solve_context);
            break;
          default:
            AssertThrow(false, dealii::ExcMessage("Unknown solver type"));
        }
      solvers.push_back(std::move(solver));
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
  grid_refiner.add_refinement_marker(std::make_shared<NucleusRefinementFunction<dim>>(
    user_inputs.get_nucleation_parameters(),
    pf_tools->nuclei_list));
  dealii::types::global_dof_index old_dofs = dof_manager.get_total_dofs();
  dealii::types::global_dof_index new_dofs = 0;
  for (unsigned int remesh_index = 0;
       remesh_index < (user_inputs.get_spatial_discretization().max_refinement -
                       user_inputs.get_spatial_discretization().min_refinement);
       remesh_index++)
    {
      // Perform grid refinement
      ConditionalOStreams::pout_base() << "performing grid refinement...\n" << std::flush;
      Timer::start_section("Grid refinement");
      grid_refiner.do_adaptive_refinement(solvers);
      Timer::end_section("Grid refinement");

      // Update the ghosts
      Timer::start_section("Update ghosts");
      for (auto &solver : solvers)
        {
          solver->update_ghosts();
        }
      Timer::end_section("Update ghosts");

      // Recalculate the total DoFs
      new_dofs = dof_manager.get_total_dofs();

      // Check for convergence
      if (old_dofs == new_dofs)
        {
          break;
        }
      old_dofs = new_dofs;
    }

  // InvM
  solve_context.get_invm_manager().compute_invm();
}

template <unsigned int dim, unsigned int degree, typename number>
void
Problem<dim, degree, number>::solve()
{
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

  ConditionalOStreams::pout_base() << "\n";

  ConditionalOStreams::pout_summary()
    << "================================================\n"
       "  Solve\n"
    << "================================================\n"
    << std::flush;

  const UserInputParameters<dim> &user_inputs = *user_inputs_ptr;
  const TemporalDiscretization   &time_info   = user_inputs.get_temporal_discretization();
  SimulationTimer                &sim_timer   = solve_context.get_simulation_timer();
  // Main time-stepping loop
  while (sim_timer.get_increment() < time_info.num_increments)
    {
      // Solve a single increment
      // Includes nucleation, refinement, constraints, solve, output, and update
      solve_increment(sim_timer);
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

  // Print timer summary
  Timer::print_summary();
}

template <unsigned int dim, unsigned int degree, typename number>
void
Problem<dim, degree, number>::solve_increment(SimulationTimer &sim_timer)
{
  const UserInputParameters<dim> &user_inputs = *user_inputs_ptr;
  unsigned int                    increment   = sim_timer.get_increment();
  bool is_output_increment = user_inputs.get_output_parameters().should_output(increment);
  bool is_nucleation_increment =
    user_inputs.get_nucleation_parameters().should_attempt_nucleation(increment);

  // Check for stochastic nucleation
  Timer::start_section("Check for nucleation");
  bool any_nucleation_occurred =
    is_nucleation_increment &&
    NucleationManager<dim, degree, number>::attempt_nucleation(solve_context,
                                                               pf_tools->nuclei_list);

  Timer::end_section("Check for nucleation");

  // Perform grid refinement if necessary
  if (user_inputs.get_spatial_discretization().has_adaptivity &&
      (user_inputs.get_spatial_discretization().should_refine_mesh(increment) ||
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

  // Update the time-dependent constraints
  Timer::start_section("Update time-dependent constraints");
  // TODO: Loop over levels, pass in current time
  constraint_manager.update_time_dependent_constraints(field_attributes);
  Timer::end_section("Update time-dependent constraints");

  // Solve a single increment
  Timer::start_section("Solve Increment");
  for (auto &solver : solvers)
    {
      SolveTiming solve_timing = solver->get_solve_group().solve_timing;
      if ((solve_timing == PostProcess && !is_output_increment) ||
          (solve_timing == NucleationRate &&
           !(is_nucleation_increment || is_output_increment)))
        {
          continue;
        }
      solver->solve();
    }
  Timer::end_section("Solve Increment");

  // Output results if needed
  if (user_inputs.get_output_parameters().should_output(increment))
    {
      Timer::start_section("Output");
      SolutionOutput<dim, number>(field_attributes,
                                  solve_context.get_solution_indexer(),
                                  sim_timer,
                                  dof_manager,
                                  degree,
                                  "solution",
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
                     dof_manager.get_dof_handler(index),
                     solution)
                << "\n";
            }
          else if (tensor_rank == TensorRank::Scalar)
            {
              ConditionalOStreams::pout_base()
                << Integrator<dim, degree, number>::template integrate<0>(
                     dof_manager.get_dof_handler(index),
                     solution)
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
}

#include "core/problem.inst"

PRISMS_PF_END_NAMESPACE