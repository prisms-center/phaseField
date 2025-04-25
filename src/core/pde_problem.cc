// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>
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

#include <prismspf/solvers/explicit_constant_solver.h>
#include <prismspf/solvers/explicit_postprocess_solver.h>
#include <prismspf/solvers/explicit_solver.h>
#include <prismspf/solvers/nonexplicit_auxiliary_solver.h>
#include <prismspf/solvers/nonexplicit_linear_solver.h>
#include <prismspf/solvers/nonexplicit_self_nonlinear_solver.h>

#include <prismspf/utilities/compute_integral.h>
#include <prismspf/utilities/element_volume.h>

#include <prismspf/config.h>

#include <memory>
#include <mpi.h>
#include <ostream>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
PDEProblem<dim, degree>::PDEProblem(
  const userInputParameters<dim>                                &_user_inputs,
  const std::shared_ptr<const PDEOperator<dim, degree, double>> &_pde_operator,
  const std::shared_ptr<const PDEOperator<dim, degree, float>>  &_pde_operator_float)
  : user_inputs(&_user_inputs)
  , mg_info(_user_inputs)
  , triangulation_handler(_user_inputs, mg_info)
  , constraint_handler(_user_inputs, mg_info)
  , matrix_free_handler()
  , multigrid_matrix_free_handler(0, 0)
  , invm_handler(*_user_inputs.var_attributes)
  , solution_handler(*_user_inputs.var_attributes, mg_info)
  , dof_handler(_user_inputs, mg_info)
  , explicit_constant_solver(_user_inputs,
                             matrix_free_handler,
                             invm_handler,
                             constraint_handler,
                             dof_handler,
                             mapping,
                             solution_handler,
                             _pde_operator)
  , explicit_solver(_user_inputs,
                    matrix_free_handler,
                    invm_handler,
                    constraint_handler,
                    dof_handler,
                    mapping,
                    solution_handler,
                    _pde_operator)
  , postprocess_explicit_solver(_user_inputs,
                                matrix_free_handler,
                                invm_handler,
                                constraint_handler,
                                dof_handler,
                                mapping,
                                solution_handler,
                                _pde_operator)
  , nonexplicit_auxiliary_solver(_user_inputs,
                                 matrix_free_handler,
                                 triangulation_handler,
                                 invm_handler,
                                 constraint_handler,
                                 dof_handler,
                                 mapping,
                                 multigrid_matrix_free_handler,
                                 solution_handler,
                                 _pde_operator)
  , nonexplicit_linear_solver(_user_inputs,
                              matrix_free_handler,
                              triangulation_handler,
                              invm_handler,
                              constraint_handler,
                              dof_handler,
                              mapping,
                              multigrid_matrix_free_handler,
                              solution_handler,
                              _pde_operator,
                              _pde_operator_float,
                              mg_info)
  , nonexplicit_self_nonlinear_solver(_user_inputs,
                                      matrix_free_handler,
                                      triangulation_handler,
                                      invm_handler,
                                      constraint_handler,
                                      dof_handler,
                                      mapping,
                                      multigrid_matrix_free_handler,
                                      solution_handler,
                                      _pde_operator,
                                      _pde_operator_float,
                                      mg_info)
{}

template <unsigned int dim, unsigned int degree>
void
PDEProblem<dim, degree>::init_system()
{
  timer::serial_timer().enter_subsection("Initialization");

  const unsigned int n_proc = dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
  conditionalOStreams::pout_base() << "number of processes: " << n_proc << "\n"
                                   << std::flush;

  const unsigned int n_vect_doubles = dealii::VectorizedArray<double>::size();
  const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;

  conditionalOStreams::pout_base()
    << "vectorization over " << n_vect_doubles << " doubles = " << n_vect_bits
    << " bits (" << dealii::Utilities::System::get_current_vectorization_level() << ')'
    << "\n"
    << std::flush;

  // Create the SCALAR/VECTOR FESystem's, if applicable
  conditionalOStreams::pout_base() << "creating FESystem...\n" << std::flush;
  for (const auto &[index, variable] : *user_inputs->var_attributes)
    {
      if (variable.field_type == fieldType::SCALAR &&
          fe_system.find(fieldType::SCALAR) == fe_system.end())
        {
          fe_system.emplace(fieldType::SCALAR,
                            dealii::FESystem<dim>(dealii::FE_Q<dim>(
                                                    dealii::QGaussLobatto<1>(degree + 1)),
                                                  1));
          conditionalOStreams::pout_summary() << "  made FESystem for scalar fields\n"
                                              << std::flush;
        }
      else if (variable.field_type == fieldType::VECTOR &&
               fe_system.find(fieldType::VECTOR) == fe_system.end())
        {
          fe_system.emplace(fieldType::VECTOR,
                            dealii::FESystem<dim>(dealii::FE_Q<dim>(
                                                    dealii::QGaussLobatto<1>(degree + 1)),
                                                  dim));
          conditionalOStreams::pout_summary() << "  made FESystem for vector fields\n"
                                              << std::flush;
        }
    }

  // Create the mesh
  conditionalOStreams::pout_base() << "creating triangulation...\n" << std::flush;
  triangulation_handler.generate_mesh();

  // Print multigrid info
  mg_info.print();

  // Create the dof handlers.
  conditionalOStreams::pout_base() << "creating DoFHandlers...\n" << std::flush;
  dof_handler.init(triangulation_handler, fe_system, mg_info);

  // Create the constraints
  conditionalOStreams::pout_base() << "creating constraints...\n" << std::flush;
  constraint_handler.make_constraints(mapping, dof_handler.get_dof_handlers());
  if (mg_info.has_multigrid())
    {
      const unsigned int min_level = mg_info.get_mg_min_level();
      const unsigned int max_level = mg_info.get_mg_max_level();
      for (unsigned int level = min_level; level <= max_level; ++level)
        {
          conditionalOStreams::pout_base()
            << "creating multigrid constraints at level " << level << "...\n"
            << std::flush;
          constraint_handler.make_mg_constraints(mapping,
                                                 dof_handler.get_mg_dof_handlers(level),
                                                 level);
        }
    }

  // Reinit the matrix-free objects
  conditionalOStreams::pout_base() << "initializing matrix-free objects...\n"
                                   << std::flush;
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
          conditionalOStreams::pout_base()
            << "initializing multgrid matrix-free object at level " << level << "...\n"
            << std::flush;
          multigrid_matrix_free_handler[level].reinit(
            mapping,
            dof_handler.get_mg_dof_handlers(level),
            constraint_handler.get_mg_constraints(level),
            dealii::QGaussLobatto<1>(degree + 1));
        }
    }

  // Initialize the solution set
  conditionalOStreams::pout_base() << "initializing solution set...\n" << std::flush;
  solution_handler.init(matrix_free_handler);
  if (mg_info.has_multigrid())
    {
      solution_handler.mg_init(multigrid_matrix_free_handler);
    }

  // Initialize the invm and compute it
  // TODO (landinjm): Output the invm for debug mode. This will create a lot of bloat in
  // the output directory so we should create a separate flag and/or directory for this.
  conditionalOStreams::pout_base() << "initializing invm...\n" << std::flush;
  invm_handler.initialize(matrix_free_handler.get_matrix_free());
  invm_handler.compute_invm();

  // Initialize the element volumes and compute them
  // TODO (landinjm): Output the element volumes for debug mode. This will create a lot of
  // bloat in the output directory so we should create a separate flag and/or directory
  // for this.
  conditionalOStreams::pout_base() << "initializing element volumes...\n" << std::flush;
  element_volume.initialize(matrix_free_handler.get_matrix_free());
  element_volume.compute_element_volume(fe_system.begin()->second);

  // Initialize the solver types
  conditionalOStreams::pout_base() << "initializing solvers...\n" << std::flush;

  explicit_constant_solver.init();
  explicit_solver.init();
  postprocess_explicit_solver.init();
  nonexplicit_auxiliary_solver.init();
  nonexplicit_linear_solver.init();
  nonexplicit_self_nonlinear_solver.init();

  // Update ghosts
  solution_handler.update_ghosts();

  // Solve the auxiliary fields at the 0th step
  conditionalOStreams::pout_base() << "solving auxiliary variables in 0th timestep...\n"
                                   << std::flush;
  nonexplicit_auxiliary_solver.solve();
  solution_handler.update_ghosts();

  // Solve the linear time-independent fields at the 0th step
  conditionalOStreams::pout_base()
    << "solving linear time-independent variables in 0th timestep...\n"
    << std::flush;
  nonexplicit_linear_solver.solve();
  solution_handler.update_ghosts();

  // Solve the self-nonlinear time-independent fields at the 0th step
  conditionalOStreams::pout_base()
    << "solving self-nonlinear time-independent variables in 0th timestep...\n"
    << std::flush;
  nonexplicit_self_nonlinear_solver.solve();
  solution_handler.update_ghosts();

  // Solve the postprocessed fields at the 0th step
  conditionalOStreams::pout_base()
    << "solving postprocessed variables in 0th timestep...\n"
    << std::flush;
  postprocess_explicit_solver.solve();
  solution_handler.update_ghosts();

  // Output initial condition
  conditionalOStreams::pout_base() << "outputting initial condition...\n" << std::flush;
  solutionOutput<dim>(solution_handler.get_solution_vector(),
                      dof_handler.get_dof_handlers(),
                      degree,
                      "solution",
                      *user_inputs);

  timer::serial_timer().leave_subsection();
}

template <unsigned int dim, unsigned int degree>
void
PDEProblem<dim, degree>::solve_increment()
{
  timer::serial_timer().enter_subsection("Solve Increment");

  // Update ghosts
  solution_handler.update_ghosts();

  // TOOD (landinjm): I think I have to update the ghosts after each solve. This should be
  // apparent in an application that includes multiple of these solve types. Also only
  // update ghosts that need to be. It's wasteful to over communicate.

  // Solve a single increment
  explicit_solver.solve();
  solution_handler.update_ghosts();

  nonexplicit_auxiliary_solver.solve();
  solution_handler.update_ghosts();

  nonexplicit_linear_solver.solve();
  solution_handler.update_ghosts();

  nonexplicit_self_nonlinear_solver.solve();
  solution_handler.update_ghosts();

  timer::serial_timer().leave_subsection();
}

template <unsigned int dim, unsigned int degree>
void
PDEProblem<dim, degree>::solve()
{
  conditionalOStreams::pout_summary()
    << "================================================\n"
       "  Initialization\n"
    << "================================================\n"
    << std::flush;

  CALI_MARK_BEGIN("Initialization");
  init_system();
  CALI_MARK_END("Initialization");

  conditionalOStreams::pout_base() << "\n";

  conditionalOStreams::pout_summary()
    << "================================================\n"
       "  Solve\n"
    << "================================================\n"
    << std::flush;
  while (user_inputs->temporal_discretization.get_current_increment() <
         user_inputs->temporal_discretization.get_total_increments())
    {
      user_inputs->temporal_discretization.update_current_increment();
      user_inputs->temporal_discretization.update_current_time();

      CALI_MARK_BEGIN("Solve Increment");
      solve_increment();
      CALI_MARK_END("Solve Increment");

      if (user_inputs->output_parameters.should_output(
            user_inputs->temporal_discretization.get_current_increment()))
        {
          // Ideally just update ghosts that need to be here.
          solution_handler.update_ghosts();
          postprocess_explicit_solver.solve();

          // TODO (landinjm): Do I need to zero out the ghost values when outputting the
          // solution?
          solutionOutput<dim>(solution_handler.get_solution_vector(),
                              dof_handler.get_dof_handlers(),
                              degree,
                              "solution",
                              *user_inputs);

          // Update the ghost again so we can call compute integral. May as well wrap this
          // into computeIntegral.
          solution_handler.update_ghosts();

          // Print the l2-norms and integrals of each solution
          conditionalOStreams::pout_base()
            << "Iteration: "
            << user_inputs->temporal_discretization.get_current_increment() << "\n";
          for (const auto &[index, vector] : solution_handler.get_solution_vector())
            {
              conditionalOStreams::pout_base()
                << "  Solution index " << index << " l2-norm: " << vector->l2_norm()
                << " integrated value: ";

              const auto local_field_type =
                user_inputs->var_attributes->at(index).field_type;

              if (local_field_type == fieldType::VECTOR)
                {
                  std::vector<double> integrated_values(dim, 0.0);
                  integral_computer.compute_integral(
                    integrated_values,
                    *dof_handler.get_dof_handlers()[index],
                    *vector);

                  for (unsigned int dimension = 0; dimension < dim; dimension++)
                    {
                      conditionalOStreams::pout_base()
                        << integrated_values[dimension] << " ";
                    }
                }
              else
                {
                  double integrated_value = 0.0;
                  integral_computer.compute_integral(
                    integrated_value,
                    *dof_handler.get_dof_handlers()[index],
                    *vector);

                  conditionalOStreams::pout_base() << integrated_value;
                }

              conditionalOStreams::pout_base() << "\n";
            }
          conditionalOStreams::pout_base() << "\n" << std::flush;
        }
    }
}

template <unsigned int dim, unsigned int degree>
void
PDEProblem<dim, degree>::run()
{
  solve();

#ifndef PRISMS_PF_WITH_CALIPER
  timer::print_summary();
#endif
}

INSTANTIATE_BI_TEMPLATE(PDEProblem)

PRISMS_PF_END_NAMESPACE
