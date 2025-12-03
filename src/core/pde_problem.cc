// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/grid_refiner.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/pde_problem.h>
#include <prismspf/core/phase_field_tools.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/solution_output.h>
#include <prismspf/core/solver_handler.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/concurrent_constant_solver.h>
#include <prismspf/solvers/concurrent_explicit_postprocess_solver.h>
#include <prismspf/solvers/concurrent_explicit_solver.h>
#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/sequential_auxiliary_solver.h>
#include <prismspf/solvers/sequential_linear_solver.h>
#include <prismspf/solvers/sequential_self_nonlinear_solver.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/utilities/element_volume.h>
#include <prismspf/utilities/integrator.h>

#include <prismspf/config.h>
#include <prismspf/nucleation/nucleation.h>
#include <prismspf/nucleation/nucleus_refinement_function.h>

#include "prismspf/user_inputs/temporal_discretization.h"

#include <memory>
#include <mpi.h>
#include <ostream>
#include <vector>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali_macros.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
PDEProblem<dim, degree, number>::PDEProblem(
  const UserInputParameters<dim>                                &_user_inputs,
  PhaseFieldTools<dim>                                          &_pf_tools,
  const std::shared_ptr<const PDEOperator<dim, degree, number>> &_pde_operator,
  const std::shared_ptr<const PDEOperator<dim, degree, float>>  &_pde_operator_float)
  : user_inputs(&_user_inputs)
  , pf_tools(&_pf_tools)
  , mg_info(_user_inputs)
  , triangulation_handler(_user_inputs, mg_info)
  , constraint_handler(_user_inputs, mg_info, _pde_operator, _pde_operator_float)
  , matrix_free_container(mg_info)
  , invm_handler(_user_inputs.get_variable_attributes())
  , solution_handler(_user_inputs.get_variable_attributes(), mg_info)
  , dof_handler(_user_inputs, mg_info)
  , element_volume_container(mg_info)
  , solver_context(_user_inputs,
                   matrix_free_container,
                   triangulation_handler,
                   invm_handler,
                   constraint_handler,
                   dof_handler,
                   mapping,
                   element_volume_container,
                   mg_info,
                   solution_handler,
                   _pde_operator,
                   _pde_operator_float)
  , grid_refiner_context(_user_inputs,
                         _pf_tools,
                         triangulation_handler,
                         constraint_handler,
                         matrix_free_container,
                         invm_handler,
                         solution_handler,
                         dof_handler,
                         fe_system,
                         mapping,
                         element_volume_container,
                         mg_info)
  , grid_refiner(grid_refiner_context)
  , solver_handler(solver_context)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
PDEProblem<dim, degree, number>::init_system()
{
  const unsigned int n_proc = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  ConditionalOStreams::pout_base() << "number of processes: " << n_proc << "\n"
                                   << std::flush;

  const unsigned int n_vect_doubles = dealii::VectorizedArray<number>::size();
  const unsigned int n_vect_bits    = 8 * sizeof(number) * n_vect_doubles;

  ConditionalOStreams::pout_base()
    << "vectorization over " << n_vect_doubles << " doubles = " << n_vect_bits
    << " bits (" << dealii::Utilities::System::get_current_vectorization_level() << ')'
    << "\n"
    << std::flush;

  // Create the Scalar/Vector FESystem's, if applicable
  ConditionalOStreams::pout_base() << "creating FESystem...\n" << std::flush;
  Timer::start_section("Create FESystem");
  for (const auto &[index, variable] : user_inputs->get_variable_attributes())
    {
      if (variable.field_info.tensor_rank == FieldInfo::TensorRank::Scalar &&
          fe_system.find(FieldInfo::TensorRank::Scalar) == fe_system.end())
        {
          fe_system.emplace(FieldInfo::TensorRank::Scalar,
                            dealii::FESystem<dim>(dealii::FE_Q<dim>(
                                                    dealii::QGaussLobatto<1>(degree + 1)),
                                                  1));
          ConditionalOStreams::pout_summary() << "  made FESystem for scalar fields\n"
                                              << std::flush;
        }
      else if (variable.field_info.tensor_rank == FieldInfo::TensorRank::Vector &&
               fe_system.find(FieldInfo::TensorRank::Vector) == fe_system.end())
        {
          fe_system.emplace(FieldInfo::TensorRank::Vector,
                            dealii::FESystem<dim>(dealii::FE_Q<dim>(
                                                    dealii::QGaussLobatto<1>(degree + 1)),
                                                  dim));
          ConditionalOStreams::pout_summary() << "  made FESystem for vector fields\n"
                                              << std::flush;
        }
    }
  Timer::end_section("Create FESystem");

  // Create the mesh
  ConditionalOStreams::pout_base() << "creating triangulation...\n" << std::flush;
  Timer::start_section("Generate mesh");
  triangulation_handler.generate_mesh();
  Timer::end_section("Generate mesh");

  // Print multigrid info
  mg_info.print();

  // Create the dof handlers.
  ConditionalOStreams::pout_base() << "creating DoFHandlers...\n" << std::flush;
  Timer::start_section("reinitialize DoFHandlers");
  dof_handler.init(triangulation_handler, fe_system, mg_info);
  Timer::end_section("reinitialize DoFHandlers");

  // Create the constraints
  ConditionalOStreams::pout_base() << "creating constraints...\n" << std::flush;
  Timer::start_section("Create constraints");
  constraint_handler.make_constraints(mapping, dof_handler.get_dof_handlers());
  if (mg_info.has_multigrid())
    {
      const unsigned int min_level = mg_info.get_mg_min_level();
      const unsigned int max_level = mg_info.get_mg_max_level();
      for (unsigned int level = min_level; level <= max_level; ++level)
        {
          ConditionalOStreams::pout_base()
            << "creating multigrid constraints at level " << level << "...\n"
            << std::flush;
          constraint_handler.make_mg_constraints(mapping,
                                                 dof_handler.get_mg_dof_handlers(level),
                                                 level);
        }
    }
  Timer::end_section("Create constraints");

  // Reinit the matrix-free objects
  matrix_free_container.template reinit<degree, 1>(mapping,
                                                   dof_handler,
                                                   constraint_handler,
                                                   dealii::QGaussLobatto<1>(degree + 1));

  // reinitialize the solution set
  solution_handler.init(matrix_free_container);

  // reinitialize the invm and compute it
  // TODO (landinjm): Output the invm for debug mode. This will create a lot of bloat in
  // the output directory so we should create a separate flag and/or directory for this.
  ConditionalOStreams::pout_base() << "initializing invm...\n" << std::flush;
  Timer::start_section("reinitialize invm");
  invm_handler.initialize(matrix_free_container.get_matrix_free());
  invm_handler.compute_invm();
  Timer::end_section("reinitialize invm");

  // reinitialize the element volumes and compute them
  element_volume_container.initialize(matrix_free_container);
  element_volume_container.compute_element_volume();

  // Initialize the solver types
  solver_handler.init();

  // Update the ghosts
  Timer::start_section("Update ghosts");
  solution_handler.update_ghosts();
  Timer::end_section("Update ghosts");

  // Perform the initial grid refinement. For this one, we have to do a loop to sufficient
  // coarsen cells to the minimum level
  ConditionalOStreams::pout_base() << "initializing grid refiner..." << std::flush;
  grid_refiner.init(fe_system);
  grid_refiner.add_refinement_marker(std::make_shared<NucleusRefinementFunction<dim>>(
    user_inputs->get_nucleation_parameters(),
    pf_tools->nuclei_list));
  dealii::types::global_dof_index old_dofs = dof_handler.get_total_dofs();
  dealii::types::global_dof_index new_dofs = 0;
  for (unsigned int remesh_index = 0;
       remesh_index < (user_inputs->get_spatial_discretization().get_max_refinement() -
                       user_inputs->get_spatial_discretization().get_min_refinement());
       remesh_index++)
    {
      // Perform grid refinement
      ConditionalOStreams::pout_base() << "performing grid refinement...\n" << std::flush;
      Timer::start_section("Grid refinement");
      grid_refiner.do_adaptive_refinement();
      Timer::end_section("Grid refinement");

      // Reinitialize the solver types
      solver_handler.reinit();

      // Update the ghosts
      Timer::start_section("Update ghosts");
      solution_handler.update_ghosts();
      Timer::end_section("Update ghosts");

      // Recalculate the total DoFs
      new_dofs = dof_handler.get_total_dofs();

      if (old_dofs == new_dofs)
        {
          break;
        }
      old_dofs = new_dofs;
    }

  // Solve the 0th timestep
  solver_handler.solve(0, true);

  // Output initial condition
  Timer::start_section("Output");
  ConditionalOStreams::pout_base() << "outputting initial condition...\n" << std::flush;
  SolutionOutput<dim, number>(solution_handler.get_solution_vector(),
                              dof_handler.get_dof_handlers(),
                              degree,
                              "solution",
                              *user_inputs);

  // Print the l2-norms and integrals of each solution
  ConditionalOStreams::pout_base()
    << "Iteration: " << user_inputs->get_temporal_discretization().get_increment()
    << "\n";
  for (const auto &[index, vector] : solution_handler.get_solution_vector())
    {
      ConditionalOStreams::pout_base()
        << "  Solution index " << index << " l2-norm: " << vector->l2_norm()
        << " integrated value: ";

      const auto local_field_type =
        user_inputs->get_variable_attributes().at(index).field_info.tensor_rank;

      if (local_field_type == FieldInfo::TensorRank::Vector)
        {
          std::vector<number> integrated_values(dim, 0.0);
          integrator.compute_integral(integrated_values,
                                      *dof_handler.get_dof_handlers()[index],
                                      *vector);

          for (unsigned int dimension = 0; dimension < dim; dimension++)
            {
              ConditionalOStreams::pout_base() << integrated_values[dimension] << " ";
            }
        }
      else
        {
          number integrated_value = 0.0;
          integrator.compute_integral(integrated_value,
                                      *dof_handler.get_dof_handlers()[index],
                                      *vector);

          ConditionalOStreams::pout_base() << integrated_value;
        }

      ConditionalOStreams::pout_base() << "\n";
    }
  ConditionalOStreams::pout_base() << "\n" << std::flush;
  Timer::end_section("Output");
}

template <unsigned int dim, unsigned int degree, typename number>
void
PDEProblem<dim, degree, number>::solve_increment()
{
  Timer::start_section("Update time-dependent constraints");
  unsigned int increment = user_inputs->get_temporal_discretization().get_increment();

  // Update the time-dependent constraints
  if (user_inputs->get_boundary_parameters().has_time_dependent_bcs())
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
    user_inputs->get_spatial_discretization().should_refine_mesh(increment) ||
    user_inputs->get_output_parameters().should_output(increment) ||
    (user_inputs->get_nucleation_parameters().postprocessed_nucleation_rate_exists() &&
     user_inputs->get_nucleation_parameters().should_attempt_nucleation(increment));

  // Solve a single increment
  solver_handler.solve(increment, update_postprocessed);
}

template <unsigned int dim, unsigned int degree, typename number>
void
PDEProblem<dim, degree, number>::solve()
{
  ConditionalOStreams::pout_summary()
    << "================================================\n"
       "  Initialization\n"
    << "================================================\n"
    << std::flush;

  Timer::start_section("Initialization");
  init_system();
  Timer::end_section("Initialization");

  ConditionalOStreams::pout_base() << "\n";

  const TemporalDiscretization &time_info = user_inputs->get_temporal_discretization();

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
  while (time_info.get_increment() < time_info.get_total_increments())
    {
      unsigned int increment               = time_info.get_increment();
      bool         any_nucleation_occurred = false;
      if (user_inputs->get_nucleation_parameters().should_attempt_nucleation(increment))
        {
          any_nucleation_occurred =
            NucleationHandler<dim, degree, number>::attempt_nucleation(
              solver_context,
              pf_tools->nuclei_list);
        }
      if (user_inputs->get_spatial_discretization().get_has_adaptivity() &&
          (user_inputs->get_spatial_discretization().should_refine_mesh(increment) ||
           any_nucleation_occurred))
        {
          // Perform grid refinement
          ConditionalOStreams::pout_base()
            << "[Increment " << time_info.get_increment() << "] : Grid Refinement\n";
          Timer::start_section("Grid refinement");
          grid_refiner.do_adaptive_refinement();
          Timer::end_section("Grid refinement");
          ConditionalOStreams::pout_base() << "\n" << std::flush;

          // Reinitialize the solver types
          solver_handler.reinit();

          // Update the ghosts
          Timer::start_section("Update ghosts");
          solution_handler.update_ghosts();
          Timer::end_section("Update ghosts");
        }

      // Update the time
      time_info.update_increment();
      time_info.update_time();
      increment = time_info.get_increment();

      Timer::start_section("Solve Increment");
      solve_increment();
      Timer::end_section("Solve Increment");

      if (user_inputs->get_output_parameters().should_output(increment))
        {
          Timer::start_section("Output");
          SolutionOutput<dim, number>(solution_handler.get_solution_vector(),
                                      dof_handler.get_dof_handlers(),
                                      degree,
                                      "solution",
                                      *user_inputs);

          // Print the l2-norms and integrals of each solution
          ConditionalOStreams::pout_base()
            << "Iteration: " << time_info.get_increment() << "\n";
          for (const auto &[index, vector] : solution_handler.get_solution_vector())
            {
              ConditionalOStreams::pout_base()
                << "  Solution index " << index << " l2-norm: " << vector->l2_norm()
                << " integrated value: ";

              const auto local_field_type =
                user_inputs->get_variable_attributes().at(index).field_info.tensor_rank;

              if (local_field_type == FieldInfo::TensorRank::Vector)
                {
                  std::vector<number> integrated_values(dim, 0.0);
                  integrator.compute_integral(integrated_values,
                                              *dof_handler.get_dof_handlers()[index],
                                              *vector);

                  for (unsigned int dimension = 0; dimension < dim; dimension++)
                    {
                      ConditionalOStreams::pout_base()
                        << integrated_values[dimension] << " ";
                    }
                }
              else
                {
                  number integrated_value = 0.0;
                  integrator.compute_integral(integrated_value,
                                              *dof_handler.get_dof_handlers()[index],
                                              *vector);

                  ConditionalOStreams::pout_base() << integrated_value;
                }

              ConditionalOStreams::pout_base() << "\n";
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
}

template <unsigned int dim, unsigned int degree, typename number>
void
PDEProblem<dim, degree, number>::run()
{
  // Print a warning if running in DEBUG mode
  ConditionalOStreams::pout_verbose()
    << "\n\n================================================\n"
       "  Warning: running in DEBUG mode \n"
    << "================================================\n\n\n"
    << std::flush;

  solve();

  Timer::print_summary();
}

#include "core/pde_problem.inst"

PRISMS_PF_END_NAMESPACE
