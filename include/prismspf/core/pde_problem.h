// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/index_map.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/solution_output.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/explicit_postprocess_solver.h>
#include <prismspf/solvers/explicit_solver.h>
#include <prismspf/solvers/nonexplicit_auxiliary_solver.h>
#include <prismspf/solvers/nonexplicit_linear_solver.h>
#include <prismspf/solvers/nonexplicit_self_nonlinear_solver.h>

#include <prismspf/utilities/element_volume.h>

#include <prismspf/config.h>

#include <map>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This is the main class that handles the construction and solving of
 * user-specified PDEs.
 */
template <int dim, int degree>
class PDEProblem
{
public:
  /**
   * \brief Constructor.
   */
  explicit PDEProblem(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Run initialization and solving steps of the given problem.
   */
  void
  run();

private:
  /**
   * \brief Main time-stepping loop that calls solve_increment, reinit_system,
   * output_results, etc...
   */
  void
  solve();

  /**
   * \brief Solve a single increment of the given PDEs.
   */
  void
  solve_increment();

  /**
   * \brief Initialize the system.
   */
  void
  init_system();

  /**
   * \brief Reinitialize the system.
   */
  void
  reinit_system();

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> &user_inputs;

  /**
   * \brief Index map.
   */
  indexMap index_map;

  /**
   * \brief Triangulation handler.
   */
  triangulationHandler<dim> triangulation_handler;

  /**
   * \brief Constraint handler.
   */
  constraintHandler<dim> constraint_handler;

  /**
   * \brief Matrix-free object handler for non-multigrid data.
   */
  matrixfreeHandler<dim> matrix_free_handler;

  /**
   * \brief Matrix-free object handler for multigrid data.
   */
  dealii::MGLevelObject<matrixfreeHandler<dim, float>> multigrid_matrix_free_handler;

  /**
   * \brief invm handler.
   */
  invmHandler<dim, degree> invm_handler;

  /**
   * \brief Solution handler.
   */
  solutionHandler<dim> solution_handler;

  /**
   * \brief DoF handler.
   */
  dofHandler<dim> dof_handler;

  /**
   * \brief Collection of finite element systems. This is just a collection of two
   * FESystem's: one for scalar fields and one for vector fields. For now they both use
   * FE_Q finite elements.
   */
  std::map<fieldType, dealii::FESystem<dim>> fe_system;

  /**
   * \brief Mappings to and from reference cell.
   */
  const dealii::MappingQ1<dim> mapping;

  /**
   * \brief Element volumes
   */
  elementVolume<dim, degree, double> element_volume;

  /**
   * \brief Explicit constant field solver class.
   */
  explicitConstantSolver<dim, degree> explicit_constant_solver;

  /**
   * \brief Explicit field solver class.
   */
  explicitSolver<dim, degree> explicit_solver;

  /**
   * \brief Postprocessed explicit field solver class.
   */
  explicitPostprocessSolver<dim, degree> postprocess_explicit_solver;

  /**
   * \brief Nonexplicit auxiliary field solver class.
   */
  nonexplicitAuxiliarySolver<dim, degree> nonexplicit_auxiliary_solver;

  /**
   * \brief Nonexplicit linear field solver class.
   */
  nonexplicitLinearSolver<dim, degree> nonexplicit_linear_solver;

  /**
   * \brief Nonexplicit self nonlinear field solver class.
   */
  nonexplicitSelfNonlinearSolver<dim, degree> nonexplicit_self_nonlinear_solver;
};

template <int dim, int degree>
PDEProblem<dim, degree>::PDEProblem(const userInputParameters<dim> &_user_inputs)
  : user_inputs(_user_inputs)
  , index_map(_user_inputs.var_attributes)
  , triangulation_handler(_user_inputs)
  , constraint_handler(_user_inputs)
  , matrix_free_handler()
  , multigrid_matrix_free_handler(0, 0)
  , invm_handler(_user_inputs.var_attributes)
  , solution_handler(_user_inputs.var_attributes)
  , dof_handler(_user_inputs)
  , explicit_constant_solver(user_inputs,
                             matrix_free_handler,
                             invm_handler,
                             constraint_handler,
                             dof_handler,
                             mapping,
                             solution_handler)
  , explicit_solver(user_inputs,
                    matrix_free_handler,
                    invm_handler,
                    constraint_handler,
                    dof_handler,
                    mapping,
                    solution_handler)
  , postprocess_explicit_solver(user_inputs,
                                matrix_free_handler,
                                invm_handler,
                                constraint_handler,
                                dof_handler,
                                mapping,
                                solution_handler)
  , nonexplicit_auxiliary_solver(user_inputs,
                                 matrix_free_handler,
                                 triangulation_handler,
                                 invm_handler,
                                 constraint_handler,
                                 dof_handler,
                                 mapping,
                                 multigrid_matrix_free_handler,
                                 solution_handler)
  , nonexplicit_linear_solver(user_inputs,
                              matrix_free_handler,
                              triangulation_handler,
                              invm_handler,
                              constraint_handler,
                              dof_handler,
                              mapping,
                              multigrid_matrix_free_handler,
                              solution_handler)
  , nonexplicit_self_nonlinear_solver(user_inputs,
                                      matrix_free_handler,
                                      triangulation_handler,
                                      invm_handler,
                                      constraint_handler,
                                      dof_handler,
                                      mapping,
                                      multigrid_matrix_free_handler,
                                      solution_handler)
{}

template <int dim, int degree>
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
  for (const auto &[index, variable] : user_inputs.var_attributes)
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

  // Create the dof handlers.
  conditionalOStreams::pout_base() << "creating DoFHandlers...\n" << std::flush;
  dof_handler.init(triangulation_handler, fe_system);

  // Create the constraints
  conditionalOStreams::pout_base() << "creating constraints...\n" << std::flush;
  constraint_handler.make_constraints(mapping, dof_handler.get_dof_handlers());
  if (triangulation_handler.has_setup_multigrid())
    {
      conditionalOStreams::pout_base() << "creating multigrid constraints...\n"
                                       << std::flush;
      constraint_handler.make_mg_constraints(mapping, dof_handler.get_mg_dof_handlers());
    }

  // Reinit the matrix-free objects
  conditionalOStreams::pout_base() << "initializing matrix-free objects...\n"
                                   << std::flush;
  matrix_free_handler.reinit(mapping,
                             dof_handler.get_dof_handlers(),
                             constraint_handler.get_constraints(),
                             dealii::QGaussLobatto<1>(degree + 1));
  if (triangulation_handler.has_setup_multigrid())
    {
      const unsigned int min_level = triangulation_handler.get_mg_min_level();
      const unsigned int max_level = triangulation_handler.get_mg_max_level();

      multigrid_matrix_free_handler.resize(min_level, max_level);

      for (unsigned int level = min_level; level <= max_level; ++level)
        {
          conditionalOStreams::pout_base()
            << "initializing multgrid matrix-free object at level " << level << "...\n"
            << std::flush;
          // multigrid_matrix_free_handler[level].reinit(mapping, );
        }
    }

  // Initialize the solution set
  conditionalOStreams::pout_base() << "initializing solution set...\n" << std::flush;
  solution_handler.init(matrix_free_handler);

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

  // Solve the linear time-independent fields at the 0th step
  conditionalOStreams::pout_base()
    << "solving linear time-independent variables in 0th timestep...\n"
    << std::flush;
  nonexplicit_linear_solver.solve();

  // Solve the self-nonlinear time-independent fields at the 0th step
  conditionalOStreams::pout_base()
    << "solving self-nonlinear time-independent variables in 0th timestep...\n"
    << std::flush;
  nonexplicit_self_nonlinear_solver.solve();

  // Solve the postprocessed fields at the 0th step
  conditionalOStreams::pout_base()
    << "solving postprocessed variables in 0th timestep...\n"
    << std::flush;
  postprocess_explicit_solver.solve();

  // Output initial condition
  conditionalOStreams::pout_base() << "outputting initial condition...\n" << std::flush;
  solutionOutput<dim>(solution_handler.get_solution_vector(),
                      dof_handler.get_dof_handlers(),
                      degree,
                      "solution",
                      user_inputs);

  timer::serial_timer().leave_subsection();
}

template <int dim, int degree>
void
PDEProblem<dim, degree>::solve_increment()
{
  timer::serial_timer().enter_subsection("Solve Increment");

  // Update ghosts
  solution_handler.update_ghosts();
  explicit_solver.solve();
  nonexplicit_auxiliary_solver.solve();
  nonexplicit_linear_solver.solve();
  nonexplicit_self_nonlinear_solver.solve();

  timer::serial_timer().leave_subsection();
}

template <int dim, int degree>
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
  while (user_inputs.temporal_discretization.increment <
         user_inputs.temporal_discretization.total_increments)
    {
      user_inputs.temporal_discretization.increment++;
      user_inputs.temporal_discretization.time += user_inputs.temporal_discretization.dt;

      CALI_MARK_BEGIN("Solve Increment");
      solve_increment();
      CALI_MARK_END("Solve Increment");

      if (user_inputs.output_parameters.should_output(
            user_inputs.temporal_discretization.increment))
        {
          postprocess_explicit_solver.solve();

          solutionOutput<dim>(solution_handler.get_solution_vector(),
                              dof_handler.get_dof_handlers(),
                              degree,
                              "solution",
                              user_inputs);

          // Print the l2-norms of each solution
          conditionalOStreams::pout_base()
            << "Iteration: " << user_inputs.temporal_discretization.increment << "\n";
          for (const auto &[index, vector] : solution_handler.get_solution_vector())
            {
              conditionalOStreams::pout_base()
                << "  Solution index " << index << " l2-norm: " << vector->l2_norm()
                << "\n";
            }
          conditionalOStreams::pout_base() << "\n" << std::flush;
        }
    }
}

template <int dim, int degree>
void
PDEProblem<dim, degree>::run()
{
  solve();

#ifndef PRISMS_PF_WITH_CALIPER
  timer::print_summary();
#endif
}

PRISMS_PF_END_NAMESPACE
