#ifndef explicit_solver_h
#define explicit_solver_h

#include <core/user_inputs/user_input_parameters.h>

/**
 * Forward declaration for user-implemented PDE class.
 */
template <int dim, int degree, typename number>
class customPDE;

/**
 * \brief This class handles the explicit solves of all explicit fields
 */
template <int dim, int degree>
class explicitSolver
{
public:
  using SystemMatrixType = customPDE<dim, degree, double>;

  /**
   * \brief Constructor.
   */
  explicitSolver(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Destructor.
   */
  ~explicitSolver();

  /**
   * \brief Initialize system.
   */
  void
  init();

  /**
   * \brief Reinitialize system.
   */
  void
  reinit();

  /**
   * \brief Solve a single update step.
   */
  void
  solve();

private:
  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> &user_inputs;

  /**
   * \brief PDE operator.
   */
  SystemMatrixType system_matrix;
};

template <int dim, int degree>
explicitSolver<dim, degree>::explicitSolver(const userInputParameters<dim> &_user_inputs)
  : user_inputs(_user_inputs)
  , system_matrix(_user_inputs)
{}

#endif