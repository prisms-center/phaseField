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
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/pde_problem.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/solution_output.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/concurrent_constant_solver.h>
#include <prismspf/solvers/concurrent_explicit_postprocess_solver.h>
#include <prismspf/solvers/concurrent_explicit_solver.h>
#include <prismspf/solvers/sequential_auxiliary_solver.h>
#include <prismspf/solvers/sequential_linear_solver.h>
#include <prismspf/solvers/sequential_self_nonlinear_solver.h>

#include <prismspf/utilities/element_volume.h>
#include <prismspf/utilities/integrator.h>

#include <prismspf/config.h>

#include <memory>
#include <mpi.h>
#include <ostream>
#include <vector>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali_macros.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
PDEProblem<dim, degree>::PDEProblem(
  const UserInputParameters<dim>                                &_user_inputs,
  const std::shared_ptr<const PDEOperator<dim, degree, double>> &_pde_operator,
  const std::shared_ptr<const PDEOperator<dim, degree, float>>  &_pde_operator_float)
  : user_inputs(&_user_inputs)
  , mg_info(_user_inputs)
  , triangulation_handler(_user_inputs, mg_info)
  , constraint_handler(_user_inputs, mg_info, _pde_operator, _pde_operator_float)
  , matrix_free_handler()
  , multigrid_matrix_free_handler(0, 0)
  , invm_handler(_user_inputs.get_variable_attributes())
  , solution_handler(_user_inputs.get_variable_attributes(), mg_info)
  , dof_handler(_user_inputs, mg_info)
  , solver_context(_user_inputs,
                   matrix_free_handler,
                   triangulation_handler,
                   invm_handler,
                   constraint_handler,
                   dof_handler,
                   mapping,
                   mg_info,
                   solution_handler,
                   multigrid_matrix_free_handler,
                   _pde_operator,
                   _pde_operator_float)
  , grid_refiner_context(_user_inputs,
                         triangulation_handler,
                         constraint_handler,
                         matrix_free_handler,
                         multigrid_matrix_free_handler,
                         invm_handler,
                         solution_handler,
                         dof_handler,
                         fe_system,
                         mapping,
                         element_volume,
                         mg_info)
  , grid_refiner(grid_refiner_context)
  , concurrent_constant_solver(solver_context)
  , concurrent_explicit_solver(solver_context)
  , concurrent_concurrent_explicit_postprocess_solver(solver_context)
  , sequential_auxiliary_solver(solver_context)
  , sequential_linear_solver(solver_context)
  , sequential_self_nonlinear_solver(solver_context)
  , sequential_co_nonlinear_solver(solver_context)
{}

template <unsigned int dim, unsigned int degree>
void
PDEProblem<dim, degree>::init_system()
{
  const unsigned int n_proc = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  ConditionalOStreams::pout_base() << "number of processes: " << n_proc << "\n"
                                   << std::flush;

  const unsigned int n_vect_doubles = dealii::VectorizedArray<double>::size();
  const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;

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
      if (variable.get_field_type() == FieldType::Scalar &&
          fe_system.find(FieldType::Scalar) == fe_system.end())
        {
          fe_system.emplace(FieldType::Scalar,
                            dealii::FESystem<dim>(dealii::FE_Q<dim>(
                                                    dealii::QGaussLobatto<1>(degree + 1)),
                                                  1));
          ConditionalOStreams::pout_summary() << "  made FESystem for scalar fields\n"
                                              << std::flush;
        }
      else if (variable.get_field_type() == FieldType::Vector &&
               fe_system.find(FieldType::Vector) == fe_system.end())
        {
          fe_system.emplace(FieldType::Vector,
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
  ConditionalOStreams::pout_base() << "initializing matrix-free objects...\n"
                                   << std::flush;
  Timer::start_section("reinitialize matrix-free objects");
  matrix_free_handler.reinit(mapping,
                             dof_handler.get_dof_handlers(),
                             constraint_handler.get_constraints(),
                             dealii::QGaussLobatto<1>(degree + 1));
  if (mg_info.has_multigrid())
    {
      const unsigned int min_level = mg_info.get_mg_min_level();
      const unsigned int max_level = mg_info.get_mg_max_level();
      multigrid_matrix_free_handler.resize(min_level, max_level);
      for (unsigned int level = min_level; level <= max_level; ++level)
        {
          ConditionalOStreams::pout_base()
            << "initializing multgrid matrix-free object at level " << level << "...\n"
            << std::flush;
          multigrid_matrix_free_handler[level].reinit(
            mapping,
            dof_handler.get_mg_dof_handlers(level),
            constraint_handler.get_mg_constraints(level),
            dealii::QGaussLobatto<1>(degree + 1));
        }
    }
  Timer::end_section("reinitialize matrix-free objects");

  // reinitialize the solution set
  ConditionalOStreams::pout_base() << "initializing solution set...\n" << std::flush;
  Timer::start_section("reinitialize solution set");
  solution_handler.init(matrix_free_handler);
  if (mg_info.has_multigrid())
    {
      solution_handler.mg_init(multigrid_matrix_free_handler);
    }
  Timer::end_section("reinitialize solution set");

  // reinitialize the invm and compute it
  // TODO (landinjm): Output the invm for debug mode. This will create a lot of bloat in
  // the output directory so we should create a separate flag and/or directory for this.
  ConditionalOStreams::pout_base() << "initializing invm...\n" << std::flush;
  Timer::start_section("reinitialize invm");
  invm_handler.initialize(matrix_free_handler.get_matrix_free());
  invm_handler.compute_invm();
  Timer::end_section("reinitialize invm");

  // reinitialize the element volumes and compute them
  // TODO (landinjm): Output the element volumes for debug mode. This will create a lot of
  // bloat in the output directory so we should create a separate flag and/or directory
  // for this.
  ConditionalOStreams::pout_base() << "initializing element volumes...\n" << std::flush;
  Timer::start_section("reinitialize element volumes");
  element_volume.initialize(matrix_free_handler.get_matrix_free());
  element_volume.compute_element_volume(fe_system.begin()->second);
  Timer::end_section("reinitialize element volumes");

  // Initialize the solver types
  ConditionalOStreams::pout_base() << "initializing solvers...\n" << std::flush;
  Timer::start_section("Solver initialization");
  ConditionalOStreams::pout_base()
    << "  trying to reinitialize concurrent constant solvers...\n"
    << std::flush;
  concurrent_constant_solver.init();
  ConditionalOStreams::pout_base()
    << "  trying to reinitialize concurrent explicit solvers...\n"
    << std::flush;
  concurrent_explicit_solver.init();
  ConditionalOStreams::pout_base()
    << "  trying to reinitialize concurrent explicit postprocess solvers...\n"
    << std::flush;
  concurrent_concurrent_explicit_postprocess_solver.init();
  ConditionalOStreams::pout_base()
    << "  trying to reinitialize sequential auxiliary solvers...\n"
    << std::flush;
  sequential_auxiliary_solver.init();
  ConditionalOStreams::pout_base()
    << "  trying to reinitialize sequential linear solvers...\n"
    << std::flush;
  sequential_linear_solver.init();
  ConditionalOStreams::pout_base()
    << "  trying to reinitialize sequential self-nonlinear solvers...\n"
    << std::flush;
  sequential_self_nonlinear_solver.init();
  ConditionalOStreams::pout_base()
    << "  trying to reinitialize sequential co-nonlinear solvers...\n"
    << std::flush;
  sequential_co_nonlinear_solver.init();
  Timer::end_section("Solver initialization");

  // Update the ghosts
  Timer::start_section("Update ghosts");
  solution_handler.update_ghosts();
  Timer::end_section("Update ghosts");

  // Perform the initial grid refinement. For this one, we have to do a loop to sufficient
  // coarsen cells to the minimum level
  ConditionalOStreams::pout_base() << "initializing grid refiner..." << std::flush;
  grid_refiner.init(fe_system);
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
      ConditionalOStreams::pout_base() << "initializing solvers...\n" << std::flush;
      Timer::start_section("Solver initialization");
      ConditionalOStreams::pout_base()
        << "  trying to reinitialize concurrent constant solvers...\n"
        << std::flush;
      concurrent_constant_solver.reinit();
      ConditionalOStreams::pout_base()
        << "  trying to reinitialize concurrent explicit solvers...\n"
        << std::flush;
      concurrent_explicit_solver.reinit();
      ConditionalOStreams::pout_base()
        << "  trying to reinitialize concurrent explicit postprocess solvers...\n"
        << std::flush;
      concurrent_concurrent_explicit_postprocess_solver.reinit();
      ConditionalOStreams::pout_base()
        << "  trying to reinitialize sequential auxiliary solvers...\n"
        << std::flush;
      sequential_auxiliary_solver.reinit();
      ConditionalOStreams::pout_base()
        << "  trying to reinitialize sequential linear solvers...\n"
        << std::flush;
      sequential_linear_solver.reinit();
      ConditionalOStreams::pout_base()
        << "  trying to reinitialize sequential self-nonlinear solvers...\n"
        << std::flush;
      sequential_self_nonlinear_solver.reinit();
      ConditionalOStreams::pout_base()
        << "  trying to reinitialize sequential co-nonlinear solvers...\n"
        << std::flush;
      sequential_co_nonlinear_solver.reinit();
      Timer::end_section("Solver initialization");

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

  // Solve the auxiliary fields at the 0th step
  ConditionalOStreams::pout_base() << "solving auxiliary variables in 0th timestep...\n"
                                   << std::flush;
  Timer::start_section("Auxiliary solver");
  sequential_auxiliary_solver.solve();
  Timer::end_section("Auxiliary solver");

  // Solve the linear time-independent fields at the 0th step
  ConditionalOStreams::pout_base()
    << "solving linear time-independent variables in 0th timestep...\n"
    << std::flush;
  Timer::start_section("Nonexplicit linear solver");
  sequential_linear_solver.solve();
  Timer::end_section("Nonexplicit linear solver");

  // Solve the self-nonlinear time-independent fields at the 0th step
  ConditionalOStreams::pout_base()
    << "solving self-nonlinear time-independent variables in 0th timestep...\n"
    << std::flush;
  Timer::start_section("Nonexplicit self-nonlinear solver");
  sequential_self_nonlinear_solver.solve();
  Timer::end_section("Nonexplicit self-nonlinear solver");

  // Solve the co-nonlinear time-independent fields at the 0th step
  ConditionalOStreams::pout_base()
    << "solving co-nonlinear time-independent variables in 0th timestep...\n"
    << std::flush;
  Timer::start_section("Nonexplicit co-nonlinear solver");
  sequential_co_nonlinear_solver.solve();
  Timer::end_section("Nonexplicit co-nonlinear solver");

  // Solve the postprocessed fields at the 0th step
  ConditionalOStreams::pout_base()
    << "solving postprocessed variables in 0th timestep...\n"
    << std::flush;
  Timer::start_section("Postprocess solver");
  concurrent_concurrent_explicit_postprocess_solver.solve();
  Timer::end_section("Postprocess solver");

  // Output initial condition
  Timer::start_section("Output");
  ConditionalOStreams::pout_base() << "outputting initial condition...\n" << std::flush;
  SolutionOutput<dim>(solution_handler.get_solution_vector(),
                      dof_handler.get_dof_handlers(),
                      degree,
                      "solution",
                      *user_inputs);

  // Print the l2-norms and integrals of each solution
  ConditionalOStreams::pout_base()
    << "Iteration: " << user_inputs->get_temporal_discretization().get_current_increment()
    << "\n";
  for (const auto &[index, vector] : solution_handler.get_solution_vector())
    {
      ConditionalOStreams::pout_base()
        << "  Solution index " << index << " l2-norm: " << vector->l2_norm()
        << " integrated value: ";

      const auto local_field_type =
        user_inputs->get_variable_attributes().at(index).get_field_type();

      if (local_field_type == FieldType::Vector)
        {
          std::vector<double> integrated_values(dim, 0.0);
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
          double integrated_value = 0.0;
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

template <unsigned int dim, unsigned int degree>
void
PDEProblem<dim, degree>::solve_increment()
{
  Timer::start_section("Update time-dependent constraints");
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

  // TOOD (landinjm): I think I have to update the ghosts after each solve. This should be
  // apparent in an application that includes multiple of these solve types. Also only
  // update ghosts that need to be. It's wasteful to over communicate.

  // Solve a single increment
  Timer::start_section("Explicit solver");
  concurrent_explicit_solver.solve();
  Timer::end_section("Explicit solver");

  Timer::start_section("Nonexplicit auxiliary solver");
  sequential_auxiliary_solver.solve();
  Timer::end_section("Nonexplicit auxiliary solver");

  Timer::start_section("Nonexplicit linear solver");
  sequential_linear_solver.solve();
  Timer::end_section("Nonexplicit linear solver");

  Timer::start_section("Nonexplicit self-nonlinear solver");
  sequential_self_nonlinear_solver.solve();
  Timer::end_section("Nonexplicit self-nonlinear solver");

  Timer::start_section("Nonexplicit co-nonlinear solver");
  sequential_co_nonlinear_solver.solve();
  Timer::end_section("Nonexplicit co-nonlinear solver");
}

template <unsigned int dim, unsigned int degree>
void
PDEProblem<dim, degree>::solve()
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

  // If the number of increments is 0, we return early
  if (user_inputs->get_temporal_discretization().get_total_increments() == 0)
    {
      return;
    }

  ConditionalOStreams::pout_summary()
    << "================================================\n"
       "  Solve\n"
    << "================================================\n"
    << std::flush;
  while (user_inputs->get_temporal_discretization().get_current_increment() <
         user_inputs->get_temporal_discretization().get_total_increments())
    {
      user_inputs->get_temporal_discretization().update_current_increment();
      user_inputs->get_temporal_discretization().update_current_time();

      Timer::start_section("Solve Increment");
      solve_increment();
      Timer::end_section("Solve Increment");

      if (user_inputs->get_spatial_discretization().should_refine_mesh(
            user_inputs->get_temporal_discretization().get_current_increment()))
        {
          // Perform grid refinement
          ConditionalOStreams::pout_base() << "performing grid refinement...\n"
                                           << std::flush;
          Timer::start_section("Grid refinement");
          grid_refiner.do_adaptive_refinement();
          Timer::end_section("Grid refinement");

          // Reinitialize the solver types
          ConditionalOStreams::pout_base() << "initializing solvers...\n" << std::flush;
          Timer::start_section("Solver initialization");
          ConditionalOStreams::pout_base()
            << "  trying to reinitialize concurrent constant solvers...\n"
            << std::flush;
          concurrent_constant_solver.reinit();
          ConditionalOStreams::pout_base()
            << "  trying to reinitialize concurrent explicit solvers...\n"
            << std::flush;
          concurrent_explicit_solver.reinit();
          ConditionalOStreams::pout_base()
            << "  trying to reinitialize concurrent explicit postprocess solvers...\n"
            << std::flush;
          concurrent_concurrent_explicit_postprocess_solver.reinit();
          ConditionalOStreams::pout_base()
            << "  trying to reinitialize sequential auxiliary solvers...\n"
            << std::flush;
          sequential_auxiliary_solver.reinit();
          ConditionalOStreams::pout_base()
            << "  trying to reinitialize sequential linear solvers...\n"
            << std::flush;
          sequential_linear_solver.reinit();
          ConditionalOStreams::pout_base()
            << "  trying to reinitialize sequential self-nonlinear solvers...\n"
            << std::flush;
          sequential_self_nonlinear_solver.reinit();
          ConditionalOStreams::pout_base()
            << "  trying to reinitialize sequential co-nonlinear solvers...\n"
            << std::flush;
          sequential_co_nonlinear_solver.reinit();
          Timer::end_section("Solver initialization");

          // Update the ghosts
          Timer::start_section("Update ghosts");
          solution_handler.update_ghosts();
          Timer::end_section("Update ghosts");
        }

      if (user_inputs->get_output_parameters().should_output(
            user_inputs->get_temporal_discretization().get_current_increment()))
        {
          Timer::start_section("Postprocess solver");
          concurrent_concurrent_explicit_postprocess_solver.solve();
          Timer::end_section("Postprocess solver");

          Timer::start_section("Output");
          SolutionOutput<dim>(solution_handler.get_solution_vector(),
                              dof_handler.get_dof_handlers(),
                              degree,
                              "solution",
                              *user_inputs);

          // Print the l2-norms and integrals of each solution
          ConditionalOStreams::pout_base()
            << "Iteration: "
            << user_inputs->get_temporal_discretization().get_current_increment() << "\n";
          for (const auto &[index, vector] : solution_handler.get_solution_vector())
            {
              ConditionalOStreams::pout_base()
                << "  Solution index " << index << " l2-norm: " << vector->l2_norm()
                << " integrated value: ";

              const auto local_field_type =
                user_inputs->get_variable_attributes().at(index).get_field_type();

              if (local_field_type == FieldType::Vector)
                {
                  std::vector<double> integrated_values(dim, 0.0);
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
                  double integrated_value = 0.0;
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
}

template <unsigned int dim, unsigned int degree>
void
PDEProblem<dim, degree>::run()
{
  solve();

  Timer::print_summary();
}

INSTANTIATE_BI_TEMPLATE(PDEProblem)

PRISMS_PF_END_NAMESPACE
