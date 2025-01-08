#ifndef pde_problem_h
#define pde_problem_h

#include <core/conditional_ostreams.h>
#include <core/matrix_free_operator.h>
#include <core/solution_output.h>
#include <core/solvers/gmg_solver.h>
#include <core/solvers/identity_solver.h>
#include <core/triangulation_handler.h>
#include <core/user_inputs/user_input_parameters.h>

/**
 * Forward declaration for user-implemented PDE class.
 */
template <int dim, int degree, typename number>
class customPDE;

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
   * \brief The field index we are solving.
   */
  const uint field_index = 0;

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> &user_inputs;

  /**
   * \brief Triangulation handler.
   */
  triangulationHandler<dim> triangulation_handler;

  /**
   * \brief The solution vector.
   */
  dealii::LinearAlgebra::distributed::Vector<double> solution;

  /**
   * \brief Triangulation DoFs.
   */
  dealii::DoFHandler<dim> dof_handler;

  /**
   * \brief Solver that we use for the linear solve
   */
  std::unique_ptr<baseSolver<dim, degree>> solver;
};

template <int dim, int degree>
PDEProblem<dim, degree>::PDEProblem(const userInputParameters<dim> &_user_inputs)
  : user_inputs(_user_inputs)
{
  if (user_inputs.linear_solve_parameters.linear_solve.at(field_index).preconditioner ==
      preconditionerType::GMG)
    {
      solver = std::make_unique<GMGSolver<dim, degree>>(user_inputs, field_index);
    }
  else
    {
      solver = std::make_unique<identitySolver<dim, degree>>(user_inputs, field_index);
    }
}

template <int dim, int degree>
void
PDEProblem<dim, degree>::init_system()
{
  solver->init_system(triangulation_handler, solution, dof_handler);
}

template <int dim, int degree>
void
PDEProblem<dim, degree>::solve_increment()
{
  solver->solve(triangulation_handler, solution, dof_handler);
}

template <int dim, int degree>
void
PDEProblem<dim, degree>::solve()
{
  triangulation_handler.generate_mesh(user_inputs);

  init_system();
  solve_increment();

  solution.update_ghost_values();
  solutionOutput<dim> output(solution, dof_handler, degree);
}

template <int dim, int degree>
void
PDEProblem<dim, degree>::run()
{
  // Print information regarding the vectorized array lanes to summary.log
  const uint n_vect_doubles = dealii::VectorizedArray<double>::size();
  const uint n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;

  conditionalOStreams::pout_summary()
    << "================================================\n"
    << "\tSolve\n"
    << "================================================\n"
    << "Vectorization over " << n_vect_doubles << " doubles = " << n_vect_bits
    << " bits (" << Utilities::System::get_current_vectorization_level() << ')' << '\n';

  // Solve solution
  solve();

  // Print summary of timing information to summary.log
}

#endif