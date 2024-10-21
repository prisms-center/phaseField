#ifndef VARIBLECONTAINER_H
#define VARIBLECONTAINER_H

#include <deal.II/base/exceptions.h>
#include <deal.II/lac/vector.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <boost/unordered_map.hpp>

#include "model_variables.h"
#include "userInputParameters.h"

/**
 * \brief This class handles the access and submission of variables and residuals.
 *
 * \tparam dim The number of dimensions in the problem.
 * \tparam degree The polynomial degree of the shape functions.
 * \tparam T The type of the scalar variable.
 */
template <int dim, int degree, typename T>
class variableContainer
{
public:
#include "typeDefs.h"
  /**
   * \brief Standard constructor for nonexplicit LHS
   */
  variableContainer(const dealii::MatrixFree<dim, double> &data,
                    const std::vector<variable_info>      &_varInfoList,
                    const std::vector<variable_info>      &_varChangeInfoList,
                    const std::vector<variable_info>      &_varOldInfoList);

  /**
   * \brief Standard constructor for explicit & nonexplicit RHS
   */
  variableContainer(const dealii::MatrixFree<dim, double> &data,
                    const std::vector<variable_info>      &_varInfoList,
                    const std::vector<variable_info>      &_varOldInfoList);

  /**
   * \brief Nonstandard constructor for postprocessing when only one index of "data"
   * should be used.
   */
  variableContainer(const dealii::MatrixFree<dim, double> &data,
                    const std::vector<variable_info>      &_varInfoList,
                    const unsigned int                    &fixed_index);

  /**
   * \brief Return the value of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   */
  T
  get_scalar_value(unsigned int global_variable_index) const;

  /**
   * \brief Return the gradient of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<1, dim, T>
  get_scalar_gradient(unsigned int global_variable_index) const;

  /**
   * \brief Return the hessian of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<2, dim, T>
  get_scalar_hessian(unsigned int global_variable_index) const;

  /**
   * \brief Return the value of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<1, dim, T>
  get_vector_value(unsigned int global_variable_index) const;

  /**
   * \brief Return the gradient of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<2, dim, T>
  get_vector_gradient(unsigned int global_variable_index) const;

  /**
   * \brief Return the hessian of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<3, dim, T>
  get_vector_hessian(unsigned int global_variable_index) const;

  /**
   * \brief Return the change in value of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   */
  T
  get_change_in_scalar_value(unsigned int global_variable_index) const;

  /**
   * \brief Return the change in gradient of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<1, dim, T>
  get_change_in_scalar_gradient(unsigned int global_variable_index) const;

  /**
   * \brief Return the change in hessian of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<2, dim, T>
  get_change_in_scalar_hessian(unsigned int global_variable_index) const;

  /**
   * \brief Return the change in value of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<1, dim, T>
  get_change_in_vector_value(unsigned int global_variable_index) const;

  /**
   * \brief Return the change in gradient of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<2, dim, T>
  get_change_in_vector_gradient(unsigned int global_variable_index) const;

  /**
   * \brief Return the change in hessian of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<3, dim, T>
  get_change_in_vector_hessian(unsigned int global_variable_index) const;

  /**
   * \brief Return the old value of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   */
  T
  get_old_scalar_value(unsigned int global_variable_index) const;

  /**
   * \brief Return the old gradient of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<1, dim, T>
  get_old_scalar_gradient(unsigned int global_variable_index) const;

  /**
   * \brief Return the old hessian of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<2, dim, T>
  get_old_scalar_hessian(unsigned int global_variable_index) const;

  /**
   * \brief Return the old value of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<1, dim, T>
  get_old_vector_value(unsigned int global_variable_index) const;

  /**
   * \brief Return the old gradient of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<2, dim, T>
  get_old_vector_gradient(unsigned int global_variable_index) const;

  /**
   * \brief Return the old hessian of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   */
  dealii::Tensor<3, dim, T>
  get_old_vector_hessian(unsigned int global_variable_index) const;

  /**
   * \brief Set the value residual of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   * \param val The value term.
   */
  void
  set_scalar_value_term_RHS(unsigned int global_variable_index, T val);

  /**
   * \brief Set the gradient residual of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   * \param grad The gradient term.
   */
  void
  set_scalar_gradient_term_RHS(unsigned int              global_variable_index,
                               dealii::Tensor<1, dim, T> grad);

  /**
   * \brief Set the value residual of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   * \param val The value term.
   */
  void
  set_vector_value_term_RHS(unsigned int              global_variable_index,
                            dealii::Tensor<1, dim, T> val);

  /**
   * \brief Set the gradient residual of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   * \param grad The gradient term.
   */
  void
  set_vector_gradient_term_RHS(unsigned int              global_variable_index,
                               dealii::Tensor<2, dim, T> grad);

  /**
   * \brief Set the value residual of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   * \param val The value term.
   */
  void
  set_scalar_value_term_LHS(unsigned int global_variable_index, T val);

  /**
   * \brief Set the gradient residual of the specified scalar field.
   *
   * \param global_variable_index The global index of the field.
   * \param grad The gradient term.
   */
  void
  set_scalar_gradient_term_LHS(unsigned int              global_variable_index,
                               dealii::Tensor<1, dim, T> grad);

  /**
   * \brief Set the value residual of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   * \param val The value term.
   */
  void
  set_vector_value_term_LHS(unsigned int              global_variable_index,
                            dealii::Tensor<1, dim, T> val);

  /**
   * \brief Set the gradient residual of the specified vector field.
   *
   * \param global_variable_index The global index of the field.
   * \param grad The gradient term.
   */
  void
  set_vector_gradient_term_LHS(unsigned int              global_variable_index,
                               dealii::Tensor<2, dim, T> grad);

  /**
   * \brief Initialize, read DOFs, and set evaulation flags for each variable.
   *
   * \param src The source vector.
   * \param cell The cell where this is done.
   */
  void
  reinit_and_eval(const std::vector<vectorType *> &src, unsigned int cell);

  /**
   * \brief Initialize, read DOFs, and set evaulation flags for the change in each
   * variable.
   *
   * \param src The source vector.
   * \param cell The cell where this is done.
   * \param var_being_solved The variable that we are evaluating.
   */
  void
  reinit_and_eval_change_in_solution(const vectorType &src,
                                     unsigned int      cell,
                                     unsigned int      var_being_solved);

  /**
   * \brief Initialize, read DOFs, and set evaulation flags for the previous timestep of
   * each variable.
   *
   * \param src The source vector.
   * \param cell The cell where this is done.
   */
  void
  reinit_and_eval_old_solution(const std::vector<vectorType *> &src, unsigned int cell);

  /**
   * Initialize the FEEvaluation object for each variable for post-processing since
   * evaluation is uneccessary.
   *
   * \param cell The cell where this is done.
   */
  void
  reinit(unsigned int cell);

  /**
   * \brief Integrate the residuals and distribute from local to global for the RHS.
   *
   * \param dst The destination vector.
   */
  void
  integrate_and_distribute(std::vector<vectorType *> &dst);

  /**
   * \brief Integrate the residuals and distribute from local to global for the LHS.
   *
   * \param dst The destination vector.
   * \param var_being_solved The variable that we are evaluating.
   */
  void
  integrate_and_distribute_change_in_solution_LHS(vectorType        &dst,
                                                  const unsigned int var_being_solved);

  /**
   * \brief The quadrature point index.
   */
  unsigned int q_point;

  /**
   * \brief Return the number of quadrature points per cell.
   */
  [[nodiscard]] unsigned int
  get_num_q_points() const;

  /**
   * \brief Return the xyz coordinates of the quadrature point.
   */
  dealii::Point<dim, T>
  get_q_point_location() const;

private:
  using scalar_FEEval = dealii::FEEvaluation<dim, degree, degree + 1, 1, double>;
  using vector_FEEval = dealii::FEEvaluation<dim, degree, degree + 1, dim, double>;

  /**
   * \brief A map of FEEvaluation objects for each active scalar variable at the current
   * timestep, specified by global variable index.
   */
  boost::unordered_map<unsigned int, std::unique_ptr<scalar_FEEval>> scalar_vars_map;

  /**
   * \brief A map of FEEvaluation objects for each active vector variable at the current
   * timestep, specified by global variable index.
   */
  boost::unordered_map<unsigned int, std::unique_ptr<vector_FEEval>> vector_vars_map;

  /**
   * \brief A map of FEEvaluation objects for the change in each active scalar variable at
   * the current timestep, specified by global variable index.
   */
  boost::unordered_map<unsigned int, std::unique_ptr<scalar_FEEval>>
    scalar_change_in_vars_map;

  /**
   * \brief A map of FEEvaluation objects for the change in each active vector variable at
   * the current timestep, specified by global variable index.
   */
  boost::unordered_map<unsigned int, std::unique_ptr<vector_FEEval>>
    vector_change_in_vars_map;

  /**
   * \brief A map of FEEvaluation objects for each active scalar variable from the
   * previous timestep, specified by global variable index.
   */
  boost::unordered_map<unsigned int, std::unique_ptr<scalar_FEEval>> scalar_old_vars_map;

  /**
   * \brief A map of FEEvaluation objects for each active vector variable from the
   * previous timestep, specified by global variable index.
   */
  boost::unordered_map<unsigned int, std::unique_ptr<vector_FEEval>> vector_old_vars_map;

  // Object containing some information about each variable

  /**
   * \brief Vector container for information about each variable (e.g., indices, whether
   * the val/grad/hess is needed, etc).
   */
  std::vector<variable_info> varInfoList;

  /**
   * \brief Vector container for information about the change in each variable (e.g.,
   * indices, whether the val/grad/hess is needed, etc).
   */
  std::vector<variable_info> varChangeInfoList;

  /**
   * \brief Vector container for information about the previous timestep of each variable
   * (e.g., indices, whether the val/grad/hess is needed, etc).
   */
  std::vector<variable_info> varOldInfoList;

  /**
   * \brief The number of fields.
   */
  unsigned int num_var;
};

#endif
