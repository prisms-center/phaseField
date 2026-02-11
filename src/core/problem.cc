// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once
#include <prismspf/core/problem.h>
#include <prismspf/core/system_wide.h>

#include "prismspf/core/simulation_timer.h"
#include "prismspf/user_inputs/temporal_discretization.h"
#include "prismspf/user_inputs/user_input_parameters.h"

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
Problem<dim, degree, number>::Problem(
  const std::vector<FieldAttributes>                            &_field_attributes,
  const std::vector<SolveGroup>                                 &_solve_groups,
  const UserInputParameters<dim>                                &_user_inputs,
  PhaseFieldTools<dim>                                          &_pf_tools,
  const std::shared_ptr<const PDEOperator<dim, degree, number>> &_pde_operator)
  : field_attributes(_field_attributes)
  , solve_groups(_solve_groups)
  , user_inputs_ptr(&_user_inputs)
  , pf_tools(&_pf_tools)
  , triangulation_manager(false)
  , dof_manager(field_attributes, solve_groups)
  , constraint_manager(field_attributes, solve_groups, dof_manager, _pde_operator)
  , solvers(solve_groups.size(), nullptr)
  , solver_context(_user_inputs,
                   triangulation_manager,
                   constraint_manager,
                   dof_manager,
                   SystemWide<dim, degree>::mapping,
                   _pde_operator)
  , grid_refiner(_user_inputs, triangulation_manager)
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
  triangulation_manager.generate_mesh();
  Timer::end_section("Generate mesh");

  // Create the dof handlers.
  ConditionalOStreams::pout_base() << "creating DoFHandlers...\n" << std::flush;
  Timer::start_section("reinitialize DoFHandlers");
  dof_manager.init(triangulation_manager, SystemWide<dim, degree>::fe_systems);
  Timer::end_section("reinitialize DoFHandlers");

  // Create the constraints
  ConditionalOStreams::pout_base() << "creating constraints...\n" << std::flush;
  Timer::start_section("Create constraints");
  for (const auto &solve_group : solve_groups)
    {
      // TODO: Loop over levels, pass in current time
      constraint_manager.make_constraints(SystemWide<dim, degree>::mapping,
                                          dof_manager.get_dof_handlers(
                                            solve_group.field_indices));
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
            solver =
              std::make_shared<ExplicitSolver<dim, degree, number>>(solve_group,
                                                                    solver_context);
            break;
          case PDEType::Linear:
            solver = std::make_shared<LinearSolver<dim, degree, number>>(solve_group,
                                                                         solver_context);
            break;
          case PDEType::Newton:
            solver = std::make_shared<NewtonSolver<dim, degree, number>>(solve_group,
                                                                         solver_context);
            break;
          default:
            AssertThrow(false, dealii::ExcMessage("Unknown solver type"));
        }
      solvers.push_back(std::move(solver));
    }
  Timer::end_section("Initialize Solvers");

  // TODO: InvM and element volume

  // Update the ghosts
  Timer::start_section("Update ghosts");
  for (auto solver : solvers)
    {
      solver.update_ghosts();
    }
  Timer::end_section("Update ghosts");

  // Perform the initial grid refinement. For this one, we have to do a loop to sufficient
  // coarsen cells to the minimum level
  ConditionalOStreams::pout_base() << "initializing grid refiner..." << std::flush;
  grid_refiner.init(SystemWide<dim, degree>::fe_systems);
  grid_refiner.add_refinement_marker(std::make_shared<NucleusRefinementFunction<dim>>(
    user_inputs.get_nucleation_parameters(),
    pf_tools->nuclei_list));
  dealii::types::global_dof_index old_dofs = dof_manager.get_total_dofs();
  dealii::types::global_dof_index new_dofs = 0;
  for (unsigned int remesh_index = 0;
       remesh_index < (user_inputs.get_spatial_discretization().get_max_refinement() -
                       user_inputs.get_spatial_discretization().get_min_refinement());
       remesh_index++)
    {
      // Perform grid refinement
      ConditionalOStreams::pout_base() << "performing grid refinement...\n" << std::flush;
      Timer::start_section("Grid refinement");
      grid_refiner.do_adaptive_refinement();
      Timer::end_section("Grid refinement");

      // Reinitialize the solvers
      for (auto solver : solvers)
        {
          solver.reinit();
        }
      // Update the ghosts
      Timer::start_section("Update ghosts");
      for (auto solver : solvers)
        {
          solver.update_ghosts();
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
  SimulationTimer                 sim_timer(time_info.get_timestep());
  // Main time-stepping loop
  while (sim_timer.get_increment() < time_info.get_total_increments())
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
  // Check for stochastic nucleation
  Timer::start_section("Check for nucleation");
  bool any_nucleation_occurred = false;
  if (user_inputs.get_nucleation_parameters().should_attempt_nucleation(increment))
    {
      any_nucleation_occurred =
        NucleationHandler<dim, degree, number>::attempt_nucleation(solver_context,
                                                                   pf_tools->nuclei_list);
    }
  Timer::end_section("Check for nucleation");

  // Perform grid refinement if necessary
  if (user_inputs.get_spatial_discretization().get_has_adaptivity() &&
      (user_inputs.get_spatial_discretization().should_refine_mesh(increment) ||
       any_nucleation_occurred))
    {
      // Perform grid refinement
      ConditionalOStreams::pout_base()
        << "[Increment " << sim_timer.get_increment() << "] : Grid Refinement\n";
      Timer::start_section("Grid refinement");
      grid_refiner.do_adaptive_refinement();
      Timer::end_section("Grid refinement");
      ConditionalOStreams::pout_base() << "\n" << std::flush;

      // Reinitialize the solvers
      for (auto solver : solvers)
        {
          solver.reinit();
        }

      // Update the ghosts
      Timer::start_section("Update ghosts");
      for (auto solver : solvers)
        {
          solver.update_ghosts();
        }
      Timer::end_section("Update ghosts");
    }

  // Update the time-dependent constraints
  Timer::start_section("Update time-dependent constraints");
  for (const auto &solve_group : solve_groups)
    {
      // TODO: Loop over levels, pass in current time
      constraint_manager.update_time_dependent_constraints(
        SystemWide<dim, degree>::mapping,
        dof_manager.get_dof_handlers(solve_group.field_indices));
    }
  Timer::end_section("Update time-dependent constraints");

  // Solve a single increment
  Timer::start_section("Solve Increment");
  for (auto solver : solvers)
    {
      solver.reinit(sim_timer);
      solver.solve();
    }

  Timer::end_section("Solve Increment");

  // Output results if needed
  if (user_inputs.get_output_parameters().should_output(increment))
    {
      Timer::start_section("Output");
      SolutionOutput<dim, number>(field_attributes,
                                  solver_context->solution_indexer,
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
            solver_context->solution_indexer.get_solution_vector(index);
          ConditionalOStreams::pout_base()
            << " Solution index " << index << " l2-norm: " << solution.l2_norm()
            << " integrated value: ";

          const auto tensor_rank = field_attributes[index].field_type;

          if (tensor_rank == FieldInfo::TensorRank::Vector)
            {
              ConditionalOStreams::pout_base()
                << Integrator<dim, degree, number>::template integrate<
                     FieldInfo::TensorRank::Vector>(dof_manager.get_dof_handler(index),
                                                    solution)
                << "\n";
            }
          else if (tensor_rank == FieldInfo::TensorRank::Scalar)
            {
              ConditionalOStreams::pout_base()
                << Integrator<dim, degree, number>::template integrate<
                     FieldInfo::TensorRank::Scalar>(dof_manager.get_dof_handler(index),
                                                    solution)
                << "\n";
            }
        }
      ConditionalOStreams::pout_base() << "\n" << std::flush;
      Timer::end_section("Output");
    }

  // Update the field labels in preparation for next increment (c_n -> c_n-1)
  for (auto solver : solvers)
    {
      solver.update();
    }
}

PRISMS_PF_END_NAMESPACE