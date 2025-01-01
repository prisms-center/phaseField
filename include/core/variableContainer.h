#ifndef VARIBLECONTAINER_H
#define VARIBLECONTAINER_H

#include <deal.II/base/exceptions.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <boost/unordered_map.hpp>

#include <core/model_variables.h>

template <int dim, int degree, typename number>
class variableContainer
{
public:
  using scalar_value = dealii::VectorizedArray<number>;
  using scalar_grad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using scalar_hess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;

  using vector_value = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using vector_grad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using vector_hess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;

  using scalar_FEEval = dealii::FEEvaluation<dim, degree, degree + 1, 1, number>;
  using vector_FEEval = dealii::FEEvaluation<dim, degree, degree + 1, dim, number>;

  // Nonexplicit LHS constructor
  variableContainer(const dealii::MatrixFree<dim, number> &data,
                    const std::vector<variable_info>      &_variable_info_list,
                    const std::vector<variable_info>      &_old_variable_info_list,
                    const std::vector<variable_info>      &_change_variable_info_list);

  // Explicit & Nonexplicit RHS constructor
  variableContainer(const dealii::MatrixFree<dim, number> &data,
                    const std::vector<variable_info>      &_variable_info_list,
                    const std::vector<variable_info>      &_old_variable_info_list);

  // Postprocess constructor
  variableContainer(const dealii::MatrixFree<dim, number> &data,
                    const std::vector<variable_info>      &_variable_info_list,
                    const unsigned int                    &fixed_index);

  [[nodiscard]] scalar_value
  get_scalar_value(unsigned int global_variable_index) const;

  [[nodiscard]] scalar_grad
  get_scalar_gradient(unsigned int global_variable_index) const;

  [[nodiscard]] scalar_hess
  get_scalar_hessian(unsigned int global_variable_index) const;

  [[nodiscard]] vector_value
  get_vector_value(unsigned int global_variable_index) const;

  [[nodiscard]] vector_grad
  get_vector_gradient(unsigned int global_variable_index) const;

  [[nodiscard]] vector_hess
  get_vector_hessian(unsigned int global_variable_index) const;

  [[nodiscard]] scalar_value
  get_change_in_scalar_value(unsigned int global_variable_index) const;

  [[nodiscard]] scalar_grad
  get_change_in_scalar_gradient(unsigned int global_variable_index) const;

  [[nodiscard]] scalar_hess
  get_change_in_scalar_hessian(unsigned int global_variable_index) const;

  [[nodiscard]] vector_value
  get_change_in_vector_value(unsigned int global_variable_index) const;

  [[nodiscard]] vector_grad
  get_change_in_vector_gradient(unsigned int global_variable_index) const;

  [[nodiscard]] vector_hess
  get_change_in_vector_hessian(unsigned int global_variable_index) const;

  [[nodiscard]] scalar_value
  get_old_scalar_value(unsigned int global_variable_index) const;

  [[nodiscard]] scalar_grad
  get_old_scalar_gradient(unsigned int global_variable_index) const;

  [[nodiscard]] scalar_hess
  get_old_scalar_hessian(unsigned int global_variable_index) const;

  [[nodiscard]] vector_value
  get_old_vector_value(unsigned int global_variable_index) const;

  [[nodiscard]] vector_grad
  get_old_vector_gradient(unsigned int global_variable_index) const;

  [[nodiscard]] vector_hess
  get_old_vector_hessian(unsigned int global_variable_index) const;

  void
  set_scalar_value_term_RHS(unsigned int global_variable_index, scalar_value val);

  void
  set_scalar_gradient_term_RHS(unsigned int global_variable_index, scalar_grad grad);

  void
  set_vector_value_term_RHS(unsigned int global_variable_index, vector_value val);

  void
  set_vector_gradient_term_RHS(unsigned int global_variable_index, vector_grad grad);

  void
  set_scalar_value_term_LHS(unsigned int global_variable_index, scalar_value val);

  void
  set_scalar_gradient_term_LHS(unsigned int global_variable_index, scalar_grad);

  void
  set_vector_value_term_LHS(unsigned int global_variable_index, vector_value val);

  void
  set_vector_gradient_term_LHS(unsigned int global_variable_index, vector_grad grad);

  // Initialize, read DOFs, and set evaulation flags for each variable
  void
  reinit_and_eval(
    const std::vector<dealii::LinearAlgebra::distributed::Vector<number> *> &src,
    unsigned int                                                             cell);

  void
  reinit_and_eval_change_in_solution(
    const dealii::LinearAlgebra::distributed::Vector<number> &src,
    unsigned int                                              cell,
    unsigned int                                              var_being_solved);

  void
  reinit_and_eval_old_solution(
    const std::vector<dealii::LinearAlgebra::distributed::Vector<number> *> &src,
    unsigned int                                                             cell);

  // Only initialize the FEEvaluation object for each variable (used for
  // post-processing)
  void
  reinit(unsigned int cell);

  // Integrate the residuals and distribute from local to global
  void
  integrate_and_distribute(
    std::vector<dealii::LinearAlgebra::distributed::Vector<number> *> &dst);

  void
  integrate_and_distribute_change_in_solution_LHS(
    dealii::LinearAlgebra::distributed::Vector<number> &dst,
    const unsigned int                                  var_being_solved);

  // The quadrature point index, a method to get the number of quadrature points
  // per cell, and a method to get the xyz coordinates for the quadrature point
  unsigned int q_point;

  [[nodiscard]] unsigned int
  get_num_q_points() const;

  [[nodiscard]] dealii::Point<dim, dealii::VectorizedArray<number>>
  get_q_point_location() const;

private:
  // Vectors of the actual FEEvaluation objects for each active variable, split
  // into scalar variables and vector variables for type reasons
  boost::unordered_map<unsigned int, std::unique_ptr<scalar_FEEval>> scalar_vars_map;
  boost::unordered_map<unsigned int, std::unique_ptr<vector_FEEval>> vector_vars_map;

  boost::unordered_map<unsigned int, std::unique_ptr<scalar_FEEval>>
    scalar_change_in_vars_map;
  boost::unordered_map<unsigned int, std::unique_ptr<vector_FEEval>>
    vector_change_in_vars_map;

  boost::unordered_map<unsigned int, std::unique_ptr<scalar_FEEval>> scalar_old_vars_map;
  boost::unordered_map<unsigned int, std::unique_ptr<vector_FEEval>> vector_old_vars_map;

  // Object containing some information about each variable (indices, whether
  // the val/grad/hess is needed, etc)
  std::vector<variable_info> variable_info_list;
  std::vector<variable_info> old_variable_info_list;
  std::vector<variable_info> change_variable_info_list;

  // The number of variables
  unsigned int num_var;
};

template class variableContainer<2, 1, double>;
template class variableContainer<3, 1, double>;

template class variableContainer<2, 2, double>;
template class variableContainer<3, 2, double>;

template class variableContainer<2, 3, double>;
template class variableContainer<3, 3, double>;

template class variableContainer<2, 4, double>;
template class variableContainer<3, 4, double>;

template class variableContainer<2, 5, double>;
template class variableContainer<3, 5, double>;

template class variableContainer<2, 6, double>;
template class variableContainer<3, 6, double>;

#endif
