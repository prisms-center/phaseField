#ifndef base_solver_h
#define base_solver_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/solver_control.h>

#include <core/boundary_conditions/constraint_handler.h>
#include <core/triangulation_handler.h>
#include <core/user_inputs/user_input_parameters.h>

/**
 * Forward declaration for user-implemented PDE class.
 */
template <int dim, int degree, typename number>
class customPDE;

/**
 * \brief Class that handles the assembly and solving of a field with the identity
 * preconditioner (no preconditioner)
 */
template <int dim, int degree>
class baseSolver
{
public:
  using SystemMatrixType = customPDE<dim, degree, double>;

  /**
   * \brief Constructor.
   */
  baseSolver(const userInputParameters<dim> &_user_inputs, const uint &_field_index);

  /**
   * \brief Destructor.
   */
  ~baseSolver() = default;

  /**
   * \brief Initialize the system.
   */
  virtual void
  init_system(const triangulationHandler<dim>                    &triangulation_handler,
              dealii::LinearAlgebra::distributed::Vector<double> &solution,
              dealii::DoFHandler<dim>                            &dof_handler) = 0;

  /**
   * \brief Reinitialize the system.
   */
  virtual void
  reinit_system() = 0;

  /**
   * \brief Solve the system Ax=b.
   */
  virtual void
  solve(const triangulationHandler<dim>                    &triangulation_handler,
        dealii::LinearAlgebra::distributed::Vector<double> &solution,
        const dealii::DoFHandler<dim>                      &dof_handler) = 0;

protected:
  /**
   * \brief Compute the solver tolerance based on the specified tolerance type.
   */
  void
  compute_solver_tolerance();

  /**
   * \brief The field index we are solving.
   */
  const uint field_index;

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> &user_inputs;

  /**
   * \brief The change in solution that is solved for to update the solution vector. This
   * is the x in Ax=b.
   */
  dealii::LinearAlgebra::distributed::Vector<double> change_in_solution;

  /**
   * \brief The residual vector. This is the b in Ax=b.
   */
  dealii::LinearAlgebra::distributed::Vector<double> residual;

  /**
   * \brief PDE operator.
   */
  SystemMatrixType system_matrix;

  /**
   * \brief Solver control.
   */
  dealii::SolverControl solver_control;

  /**
   * \brief Solver tolerance
   */
  double tolerance = 0.0;
};

template <int dim, int degree>
baseSolver<dim, degree>::baseSolver(const userInputParameters<dim> &_user_inputs,
                                    const uint                     &_field_index)
  : field_index(_field_index)
  , user_inputs(_user_inputs)
  , system_matrix(_user_inputs)
  , solver_control(
      _user_inputs.linear_solve_parameters.linear_solve.at(_field_index).max_iterations)
{}

template <int dim, int degree>
inline void
baseSolver<dim, degree>::compute_solver_tolerance()
{
  tolerance =
    user_inputs.linear_solve_parameters.linear_solve.at(field_index).tolerance_type ==
        solverToleranceType::RELATIVE_RESIDUAL_CHANGE
      ? user_inputs.linear_solve_parameters.linear_solve.at(field_index).tolerance *
          residual.l2_norm()
      : user_inputs.linear_solve_parameters.linear_solve.at(field_index).tolerance;
}

#endif