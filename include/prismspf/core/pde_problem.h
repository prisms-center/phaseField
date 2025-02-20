// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef pde_problem_h
#define pde_problem_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/element_volume.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/solution_output.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/solvers/explicit_postprocess_solver.h>
#include <prismspf/solvers/explicit_solver.h>
#include <prismspf/solvers/nonexplicit_auxiliary_solver.h>
#include <prismspf/solvers/nonexplicit_linear_solver.h>
#include <prismspf/solvers/nonexplicit_self_nonlinear_solver.h>
#include <prismspf/user_inputs/user_input_parameters.h>

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
  PDEProblem(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Destructor.
   */
  ~PDEProblem() = default;

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
  , triangulation_handler(_user_inputs)
  , constraint_handler(_user_inputs)
  , matrix_free_handler(_user_inputs)
  , multigrid_matrix_free_handler(0, 0, _user_inputs)
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

  const unsigned int n_vect_doubles = dealii::VectorizedArray<double>::size();
  const unsigned int n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;

  conditionalOStreams::pout_summary()
    << "Vectorization over " << n_vect_doubles << " doubles = " << n_vect_bits
    << " bits (" << dealii::Utilities::System::get_current_vectorization_level() << ')'
    << "\n"
    << std::flush;

  // Create the SCALAR/VECTOR FESystem's, if applicable
  conditionalOStreams::pout_summary() << "creating FESystem...\n" << std::flush;
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
  conditionalOStreams::pout_summary() << "creating triangulation...\n" << std::flush;
  triangulation_handler.generate_mesh();

  // Create the dof handlers.
  conditionalOStreams::pout_summary() << "creating DoFHandlers...\n" << std::flush;
  dof_handler.init(triangulation_handler, fe_system);

  // Create the constraints
  conditionalOStreams::pout_summary() << "creating constraints...\n" << std::flush;
  constraint_handler.make_constraints(mapping, dof_handler.dof_handlers);

  // Reinit the matrix-free objects
  conditionalOStreams::pout_summary() << "initializing matrix-free objects...\n"
                                      << std::flush;
  matrix_free_handler.reinit(mapping,
                             dof_handler.const_dof_handlers,
                             constraint_handler.get_constraints(),
                             dealii::QGaussLobatto<1>(degree + 1));

  // Initialize the solution set
  conditionalOStreams::pout_summary() << "initializing solution set...\n" << std::flush;
  solution_handler.init(matrix_free_handler);

  // Initialize the invm and compute it
  // TODO: Output the invm for debug mode
  conditionalOStreams::pout_summary() << "initializing invm...\n" << std::flush;
  invm_handler.initialize(matrix_free_handler.get_matrix_free());
  invm_handler.compute_invm();

  // Initialize the element volumes and compute them
  // TODO: Output the element volumes for debug mode
  conditionalOStreams::pout_summary() << "initializing element volumes...\n"
                                      << std::flush;
  element_volume.initialize(matrix_free_handler.get_matrix_free());
  element_volume.compute_element_volume(fe_system.begin()->second);

  // Initialize the solver types
  conditionalOStreams::pout_summary() << "initializing solvers...\n" << std::flush;

  CALI_MARK_BEGIN("Constant init");
  explicit_constant_solver.init();
  CALI_MARK_END("Constant init");

  CALI_MARK_BEGIN("Explicit init");
  explicit_solver.init();
  CALI_MARK_END("Explicit init");

  CALI_MARK_BEGIN("Postprocess init");
  postprocess_explicit_solver.init();
  CALI_MARK_END("Postprocess init");

  CALI_MARK_BEGIN("Auxiliary init");
  nonexplicit_auxiliary_solver.init();
  CALI_MARK_END("Auxiliary init");

  CALI_MARK_BEGIN("Linear init");
  nonexplicit_linear_solver.init();
  CALI_MARK_END("Linear init");

  CALI_MARK_BEGIN("Self-nonlinear init");
  nonexplicit_self_nonlinear_solver.init();
  CALI_MARK_END("Self-nonlinear init");

  // Update ghosts
  CALI_MARK_BEGIN("Update ghosts");
  solution_handler.update_ghosts();
  CALI_MARK_END("Update ghosts");

  // Solve the auxiliary fields at the 0th step
  conditionalOStreams::pout_summary()
    << "solving auxiliary variables in 0th timestep...\n"
    << std::flush;
  CALI_MARK_BEGIN("Auxiliary solve");
  nonexplicit_auxiliary_solver.solve();
  CALI_MARK_END("Auxiliary solve");

  // Solve the linear time-independent fields at the 0th step
  conditionalOStreams::pout_summary()
    << "solving linear time-independent variables in 0th timestep...\n"
    << std::flush;
  CALI_MARK_BEGIN("Linear solve");
  nonexplicit_linear_solver.solve();
  CALI_MARK_END("Linear solve");

  // Solve the self-nonlinear time-independent fields at the 0th step
  conditionalOStreams::pout_summary()
    << "solving self-nonlinear time-independent variables in 0th timestep...\n"
    << std::flush;
  CALI_MARK_BEGIN("Self-nonlinear solve");
  nonexplicit_self_nonlinear_solver.solve();
  CALI_MARK_END("Self-nonlinear solve");

  // Solve the postprocessed fields at the 0th step
  conditionalOStreams::pout_summary()
    << "solving postprocessed variables in 0th timestep...\n"
    << std::flush;
  CALI_MARK_BEGIN("Postprocess solve");
  postprocess_explicit_solver.solve();
  CALI_MARK_END("Postprocess solve");

  // Output initial condition
  conditionalOStreams::pout_summary() << "outputting initial condition...\n"
                                      << std::flush;
  solutionOutput<dim> output_solution(solution_handler.solution_set,
                                      dof_handler.const_dof_handlers,
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
  CALI_MARK_BEGIN("Update ghosts");
  solution_handler.update_ghosts();
  CALI_MARK_END("Update ghosts");

  CALI_MARK_BEGIN("Explicit solve");
  explicit_solver.solve();
  CALI_MARK_END("Explicit solve");

  CALI_MARK_BEGIN("Auxiliary solve");
  nonexplicit_auxiliary_solver.solve();
  CALI_MARK_END("Auxiliary solve");

  CALI_MARK_BEGIN("Linear solve");
  nonexplicit_linear_solver.solve();
  CALI_MARK_END("Linear solve");

  CALI_MARK_BEGIN("Self-nonlinear solve");
  nonexplicit_self_nonlinear_solver.solve();
  CALI_MARK_END("Self-nonlinear solve");

  timer::serial_timer().leave_subsection();
}

template <int dim, int degree>
void
PDEProblem<dim, degree>::solve()
{
  CALI_MARK_BEGIN("Main init");
  conditionalOStreams::pout_summary()
    << "================================================\n"
       "  Initialization\n"
    << "================================================\n"
    << std::flush;
  init_system();
  conditionalOStreams::pout_summary() << "\n";
  CALI_MARK_END("Main init");

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

      CALI_MARK_BEGIN("Solve increment");
      solve_increment();
      CALI_MARK_END("Solve increment");

      if (user_inputs.output_parameters.should_output(
            user_inputs.temporal_discretization.increment))
        {
          CALI_MARK_BEGIN("Output");

          CALI_MARK_BEGIN("Postprocess solve");
          postprocess_explicit_solver.solve();
          CALI_MARK_END("Postprocess solve");

          solutionOutput<dim> output_solution(solution_handler.solution_set,
                                              dof_handler.const_dof_handlers,
                                              degree,
                                              "solution",
                                              user_inputs);

          // Print the l2-norms of each solution
          conditionalOStreams::pout_base()
            << "Iteration: " << user_inputs.temporal_discretization.increment << "\n";
          for (const auto &[pair, vector] : solution_handler.solution_set)
            {
              conditionalOStreams::pout_base()
                << "  Solution index " << pair.first << " type " << to_string(pair.second)
                << " l2-norm: " << vector->l2_norm() << "\n";
            }
          conditionalOStreams::pout_base() << "\n" << std::flush;
          CALI_MARK_END("Output");
        }
    }
}

template <int dim, int degree>
void
PDEProblem<dim, degree>::run()
{
  CALI_MARK_BEGIN("Main solve");
  solve();
  CALI_MARK_END("Main solve");

  timer::print_summary();
}

PRISMS_PF_END_NAMESPACE

#endif