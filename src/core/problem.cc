// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once
#include <prismspf/core/problem.h>

#include "prismspf/user_inputs/user_input_parameters.h"

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
Problem<dim, degree, number>::fe_systems = {
  dealii::FESystem<dim>(dealii::FE_Q<dim>(dealii::QGaussLobatto<1>(degree + 1)), 1),
  dealii::FESystem<dim>(dealii::FE_Q<dim>(dealii::QGaussLobatto<1>(degree + 1)), dim)};

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
  , triangulation_handler(_user_inputs)
  , constraint_handler(_user_inputs, _pde_operator)
  , dof_manager(_user_inputs)
  , solver_context(_user_inputs,
                   triangulation_handler,
                   constraint_handler,
                   dof_manager,
                   mapping,
                   solution_handler,
                   _pde_operator)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
Problem<dim, degree, number>::init_system()
{
  // Initialize the solvers
  solvers.reserve(solve_groups.size());
  for (const auto &solve_group : solve_groups)
    {
      std::shared_ptr<GroupSolverBase<dim, degree, number>> solver;
      switch (solve_group.solver_type)
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
  triangulation_handler.generate_mesh();
  Timer::end_section("Generate mesh");

  // Create the dof handlers.
  ConditionalOStreams::pout_base() << "creating DoFHandlers...\n" << std::flush;
  Timer::start_section("reinitialize DoFHandlers");
  dof_manager.init(triangulation_handler, fe_systems);
  Timer::end_section("reinitialize DoFHandlers");

  // Create the constraints
  ConditionalOStreams::pout_base() << "creating constraints...\n" << std::flush;
  Timer::start_section("Create constraints");
  constraint_handler.make_constraints(mapping, dof_manager.get_dof_handlers());
  const unsigned int min_level = 0;
  const unsigned int max_level = 0;
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      ConditionalOStreams::pout_base()
        << "creating multigrid constraints at level " << level << "...\n"
        << std::flush;
      constraint_handler.make_mg_constraints(mapping,
                                             dof_manager.get_mg_dof_handlers(level),
                                             level);
    }
  Timer::end_section("Create constraints");

  // TODO: InvM and element volume

  // Update the ghosts
  Timer::start_section("Update ghosts");
  solution_handler.update_ghosts();
  Timer::end_section("Update ghosts");

  // Perform the initial grid refinement. For this one, we have to do a loop to sufficient
  // coarsen cells to the minimum level
  ConditionalOStreams::pout_base() << "initializing grid refiner..." << std::flush;
  grid_refiner.init(fe_systems);
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

  const UserInputParameters<dim> &user_inputs = *user_inputs_ptr;
  const TemporalDiscretization   &time_info   = user_inputs.get_temporal_discretization();
  SimulationTimer                 sim_timer(time_info);

  // If the number of increments is 0, we return early
  if (time_info.get_total_increments() == 0)
    {
      return;
    }

  ConditionalOStreams::pout_summary()
    << "================================================\n"
       "  Solve\n"
    << "================================================\n"
    << std::flush;
  while (sim_timer.get_increment() < time_info.get_total_increments())
    {
      unsigned int increment = sim_timer.get_increment();
      // Check for stochastic nucleation
      bool any_nucleation_occurred = false;
      if (user_inputs.get_nucleation_parameters().should_attempt_nucleation(increment))
        {
          any_nucleation_occurred =
            NucleationHandler<dim, degree, number>::attempt_nucleation(
              solver_context,
              pf_tools->nuclei_list);
        }
      // Perform grid refinement if necessary
      if (user_inputs.get_spatial_discretization().get_has_adaptivity() &&
          (user_inputs.get_spatial_discretization().should_refine_mesh(increment) ||
           any_nucleation_occurred))
        {
          // Perform grid refinement
          ConditionalOStreams::pout_base()
            << "[Increment " << time_info.get_increment() << "] : Grid Refinement\n";
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

      Timer::start_section("Solve Increment");
      solve_increment(sim_timer);
      Timer::end_section("Solve Increment");

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
            << "Iteration: " << time_info.get_increment() << "\n";
          for (unsigned int index = 0; index < field_attributes.size(); ++index)
            {
              const auto &solution =
                solver_context->solution_indexer.get_solution_vector(index);
              ConditionalOStreams::pout_base()
                << " Solution index " << index << " l2-norm: " << solution.l2_norm()
                << " integrated value: ";

              const auto tensor_rank = field_attributes[index].rank;

              if (tensor_rank == FieldInfo::TensorRank::Vector)
                {
                  ConditionalOStreams::pout_base()
                    << Integrator<dim, degree, number>::template integrate<
                         FieldInfo::TensorRank::Vector>(dof_manager.get_dof_handler(
                                                          index),
                                                        solution)
                    << "\n";
                }
              else if (tensor_rank == FieldInfo::TensorRank::Scalar)
                {
                  ConditionalOStreams::pout_base()
                    << Integrator<dim, degree, number>::template integrate<
                         FieldInfo::TensorRank::Scalar>(dof_manager.get_dof_handler(
                                                          index),
                                                        solution)
                    << "\n";
                }
            }
          ConditionalOStreams::pout_base() << "\n" << std::flush;
          Timer::end_section("Output");
        }
    }
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
  Timer::print_summary();
}

template <unsigned int dim, unsigned int degree, typename number>
void
Problem<dim, degree, number>::solve_increment()
{
  const UserInputParameters<dim> &user_inputs = *user_inputs_ptr;
  Timer::start_section("Update time-dependent constraints");
  unsigned int increment = user_inputs.get_temporal_discretization().get_increment();

  // Update the time-dependent constraints
  if (user_inputs.get_boundary_parameters().has_time_dependent_bcs())
    {
      constraint_handler
        .update_time_dependent_constraints(mapping, dof_handler.get_dof_handlers());
      if (mg_info.has_multigrid())
        {
          const unsigned int min_level = mg_info.get_mg_min_level();
          const unsigned int max_level = mg_info.get_mg_max_level();
          for (unsigned int level = min_level; level <= max_level; ++level)
            {
              constraint_handler.update_time_dependent_mg_constraints(
                mapping,
                dof_handler.get_mg_dof_handlers(level),
                level);
            }
        }
    }
  Timer::end_section("Update time-dependent constraints");

  // TODO (landinjm): I think I have to update the ghosts after each solve. This should be
  // apparent in an application that includes multiple of these solve types. Also only
  // update ghosts that need to be. It's wasteful to over communicate.
  bool update_postprocessed =
    user_inputs.get_spatial_discretization().should_refine_mesh(increment) ||
    user_inputs.get_output_parameters().should_output(increment) ||
    (user_inputs.get_nucleation_parameters().postprocessed_nucleation_rate_exists() &&
     user_inputs.get_nucleation_parameters().should_attempt_nucleation(increment));

  // Solve a single increment
  solver_handler.solve(increment, update_postprocessed);
  for (auto solver : solvers)
    {
      solver.reinit(timer);
      solver.solve();
    }
  for (auto solver : solvers)
    {
      solver.update();
    }
}

PRISMS_PF_END_NAMESPACE