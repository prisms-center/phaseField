// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class permits the access of a subset of indexed fields and gives an error
 * if any non-allowed fields are requested.
 *
 * \tparam dim The number of dimensions in the problem.
 * \tparam degree The polynomial degree of the shape functions.
 * \tparam number Datatype to use for `dealii::VectorizedArray<number>`. Either
 * double or float.
 */
template <int dim, int degree, typename number>
class variableContainer
{
public:
  using VectorType = dealii::LinearAlgebra::distributed::Vector<number>;
  using value_type = number;
  using size_type  = dealii::VectorizedArray<number>;

  /**
   * \brief Constructor.
   */
  variableContainer(
    const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
    const std::map<unsigned int, variableAttributes> &_subset_attributes,
    const std::unordered_map<std::pair<unsigned int, dependencyType>,
                             unsigned int,
                             pairHash>               &_global_to_local_solution,
    const solveType                                  &_solve_type);

  /**
   * \brief Return the value of the specified scalar field.
   */
  [[nodiscard]] size_type
  get_scalar_value(unsigned int   global_variable_index,
                   dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the gradient of the specified scalar field.
   */
  [[nodiscard]] dealii::Tensor<1, dim, size_type>
  get_scalar_gradient(unsigned int   global_variable_index,
                      dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the hessian of the specified scalar field.
   */
  [[nodiscard]] dealii::Tensor<2, dim, size_type>
  get_scalar_hessian(unsigned int   global_variable_index,
                     dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the diagonal of the hessian of the specified scalar field.
   */
  [[nodiscard]] dealii::Tensor<1, dim, size_type>
  get_scalar_hessian_diagonal(
    unsigned int   global_variable_index,
    dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the laplacian of the specified scalar field.
   */
  [[nodiscard]] size_type
  get_scalar_laplacian(unsigned int   global_variable_index,
                       dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the value of the specified vector field.
   */
  [[nodiscard]] dealii::Tensor<1, dim, size_type>
  get_vector_value(unsigned int   global_variable_index,
                   dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the gradient of the specified vector field.
   */
  [[nodiscard]] dealii::Tensor<2, dim, size_type>
  get_vector_gradient(unsigned int   global_variable_index,
                      dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the hessian of the specified vector field.
   */
  [[nodiscard]] dealii::Tensor<3, dim, size_type>
  get_vector_hessian(unsigned int   global_variable_index,
                     dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the diagonal of the hessian of the specified vector field.
   */
  [[nodiscard]] dealii::Tensor<2, dim, size_type>
  get_vector_hessian_diagonal(
    unsigned int   global_variable_index,
    dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the laplacian of the specified vector field.
   */
  [[nodiscard]] dealii::Tensor<1, dim, size_type>
  get_vector_laplacian(unsigned int   global_variable_index,
                       dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the divergence of the specified vector field.
   */
  [[nodiscard]] size_type
  get_vector_divergence(unsigned int   global_variable_index,
                        dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the symmetric gradient of the specified vector field.
   */
  [[nodiscard]] dealii::Tensor<2, dim, size_type>
  get_vector_symmetric_gradient(
    unsigned int   global_variable_index,
    dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the curl of the specified vector field. Note that this is
   * dealii::VectorizedArray<number> type for 2D and dealii::Tensor<1, dim,
   * dealii::VectorizedArray<number>> type for 3D.
   */
  [[nodiscard]] dealii::Tensor<1, (dim == 2 ? 1 : dim), size_type>
  get_vector_curl(unsigned int   global_variable_index,
                  dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Set the residual value of the specified scalar field.
   */
  void
  set_scalar_value_term(const unsigned int   &global_variable_index,
                        const size_type      &val,
                        const dependencyType &dependency_type = dependencyType::NORMAL);

  /**
   * \brief Set the residual gradient of the specified scalar field.
   */
  void
  set_scalar_gradient_term(
    const unsigned int                      &global_variable_index,
    const dealii::Tensor<1, dim, size_type> &grad,
    const dependencyType                    &dependency_type = dependencyType::NORMAL);

  /**
   * \brief Set the residual value of the specified vector field.
   */
  void
  set_vector_value_term(const unsigned int                      &global_variable_index,
                        const dealii::Tensor<1, dim, size_type> &val,
                        const dependencyType &dependency_type = dependencyType::NORMAL);

  /**
   * \brief Set the residual gradient of the specified vector field.
   */
  void
  set_vector_gradient_term(
    const unsigned int                      &global_variable_index,
    const dealii::Tensor<2, dim, size_type> &grad,
    const dependencyType                    &dependency_type = dependencyType::NORMAL);

  /**
   * \brief Apply some operator function for a given cell range and source vector to
   * some destination vector.
   */
  void
  eval_local_operator(
    const std::function<void(variableContainer &, const dealii::Point<dim, size_type> &)>
                                                &func,
    std::vector<VectorType *>                   &dst,
    const std::vector<VectorType *>             &src,
    const std::pair<unsigned int, unsigned int> &cell_range);

  /**
   * \brief Apply some operator function for a given cell range and source vector to
   * some destination vector.
   */
  void
  eval_local_operator(
    const std::function<void(variableContainer &, const dealii::Point<dim, size_type> &)>
                                                &func,
    VectorType                                  &dst,
    const std::vector<VectorType *>             &src,
    const std::pair<unsigned int, unsigned int> &cell_range);

  /**
   * \brief Apply some operator function for a given cell range and source vector to
   * some destination vector.
   */
  void
  eval_local_operator(
    const std::function<void(variableContainer &, const dealii::Point<dim, size_type> &)>
                                                &func,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::vector<VectorType *>             &src_subset,
    const std::pair<unsigned int, unsigned int> &cell_range);

  /**
   * \brief TODO (landinjm): Add comments
   */
  void
  eval_local_diagonal(
    const std::function<void(variableContainer &, const dealii::Point<dim, size_type> &)>
                                                &func,
    VectorType                                  &dst,
    const std::vector<VectorType *>             &src_subset,
    const std::pair<unsigned int, unsigned int> &cell_range);

private:
  using scalar_FEEval = dealii::FEEvaluation<dim, degree, degree + 1, 1, number>;
  using vector_FEEval = dealii::FEEvaluation<dim, degree, degree + 1, dim, number>;

  /**
   * \brief Check whether the map entry for the scalar FEEvaluation exists.
   */
  void
  scalar_FEEval_exists(const unsigned int   &dependency_index,
                       const dependencyType &dependency_type) const;

  /**
   * \brief Check whether the map entry for the vector FEEvaluation exists.
   */
  void
  vector_FEEval_exists(const unsigned int   &dependency_index,
                       const dependencyType &dependency_type) const;

  /**
   * \brief Check that a variable value/gradient/hessians was marked as needed and thus
   * properly initialized.
   */
  void
  access_valid(const unsigned int                             &dependency_index,
               const dependencyType                           &dependency_type,
               const dealii::EvaluationFlags::EvaluationFlags &flag) const;

  /**
   * \brief Check that a value is valid for submission.
   */
  void
  submission_valid(const dependencyType &dependency_type) const;

  /**
   * \brief Return the number of quadrature points.
   */
  [[nodiscard]] unsigned int
  get_n_q_points() const;

  /**
   * \brief Return the quadrate point location.
   */
  [[nodiscard]] dealii::Point<dim, size_type>
  get_q_point_location() const;

  /**
   * \brief Initialize, read DOFs, and set evaulation flags for each variable.
   */
  void
  reinit_and_eval(const std::vector<VectorType *> &src, unsigned int cell);

  /**
   * \brief Initialize, read DOFs, and set evaulation flags for each variable.
   */
  void
  reinit_and_eval(const VectorType &src, unsigned int cell);

  /**
   * \brief Initialize the cell for all dependencies of a certain variable index.
   */
  void
  reinit(unsigned int cell, const unsigned int &global_variable_index);

  /**
   * \brief Read dofs values on the cell for all dependencies of a certain variable index.
   */
  void
  read_dof_values(const std::vector<VectorType *> &src, unsigned int cell);

  /**
   * \brief Evaluate the flags on the cell for all dependencies of a certain variable
   * index.
   */
  void
  eval(const unsigned int &global_variable_index);

  /**
   * \brief Integrate the residuals and distribute from local to global.
   */
  void
  integrate_and_distribute(std::vector<VectorType *> &dst);

  /**
   * \brief Integrate the residuals and distribute from local to global.
   */
  void
  integrate_and_distribute(VectorType &dst);

  /**
   * \brief Integrate the residuals for a certain variable index.
   */
  void
  integrate(const unsigned int &global_variable_index);

  /**
   * \brief Map of FEEvaluation objects for each active scalar variables. The first
   * mapping is for the global variable and the second is for the dependencyType.
   */
  std::map<unsigned int, std::map<dependencyType, std::unique_ptr<scalar_FEEval>>>
    scalar_vars_map;

  /**
   * \brief Map of FEEvaluation objects for each active vector variables. The first
   * mapping is for the global variable and the second is for the dependencyType.
   */
  std::map<unsigned int, std::map<dependencyType, std::unique_ptr<vector_FEEval>>>
    vector_vars_map;

  /**
   * \brief The attribute list of the relevant subset of variables.
   */
  const std::map<unsigned int, variableAttributes> *subset_attributes;

  /**
   * \brief Mapping from global solution vectors to the local ones
   */
  const std::unordered_map<std::pair<unsigned int, dependencyType>,
                           unsigned int,
                           pairHash> *global_to_local_solution;

  /**
   * \brief The residual evaluation flags taken in from the subset attributes. For all
   * solve types, there is only a single unique instance of the eval flags, so we can
   * simply store it here when the constructor is called.
   */
  std::unordered_map<std::pair<unsigned int, dependencyType>,
                     dealii::EvaluationFlags::EvaluationFlags,
                     pairHash>
    src_eval_flags;

  /**
   * \brief The destination evaluation flags taken in from the subset attributes. For all
   * solve types, there is only a single unique instance of the eval flags, so we can
   * simply store it here when the constructor is called.
   */
  dealii::EvaluationFlags::EvaluationFlags dst_eval_flags =
    dealii::EvaluationFlags::EvaluationFlags::nothing;

  /**
   * \brief The solve type
   */
  solveType solve_type;

  /**
   * \brief The quadrature point index.
   */
  unsigned int q_point = 0;

  /**
   * \brief Number of DoFs per cell.
   */
  unsigned int n_dofs_per_cell = 0;

  /**
   * \brief Diagonal matrix that is used for preconditioning of scalar fields.
   */
  std::unique_ptr<dealii::AlignedVector<size_type>> scalar_diagonal;

  /**
   * \brief Diagonal matrix that is used for preconditioning of vector fields.
   */
  std::unique_ptr<dealii::AlignedVector<dealii::Tensor<1, dim, size_type>>>
    vector_diagonal;
};

PRISMS_PF_END_NAMESPACE
