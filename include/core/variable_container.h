#ifndef variable_container_h
#define variable_container_h

#include <deal.II/lac/vector.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <core/exceptions.h>
#include <core/type_enums.h>
#include <core/variable_attributes.h>

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
  using EvalFlags   = dealii::EvaluationFlags::EvaluationFlags;
  using scalarValue = dealii::VectorizedArray<number>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using scalarHess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using vectorValue = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using vectorGrad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using vectorHess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;

  /**
   * \brief Constructor.
   */
  variableContainer(const dealii::MatrixFree<dim, number> &data,
                    const AttributesList                  &_subset_attributes,
                    const solveType                       &_solve_type);

  /**
   * \brief Return the value of the specified scalar field.
   */
  scalarValue
  get_scalar_value(uint           global_variable_index,
                   dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the gradient of the specified scalar field.
   */
  scalarGrad
  get_scalar_gradient(uint           global_variable_index,
                      dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the hessian of the specified scalar field.
   */
  scalarHess
  get_scalar_hessian(uint           global_variable_index,
                     dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the diagonal of the hessian of the specified scalar field.
   */
  scalarGrad
  get_scalar_hessian_diagonal(
    uint           global_variable_index,
    dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the laplacian of the specified scalar field.
   */
  scalarValue
  get_scalar_laplacian(uint           global_variable_index,
                       dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the value of the specified vector field.
   */
  vectorValue
  get_vector_value(uint           global_variable_index,
                   dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the gradient of the specified vector field.
   */
  vectorGrad
  get_vector_gradient(uint           global_variable_index,
                      dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the hessian of the specified vector field.
   */
  vectorHess
  get_vector_hessian(uint           global_variable_index,
                     dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the diagonal of the hessian of the specified vector field.
   */
  vectorGrad
  get_vector_hessian_diagonal(
    uint           global_variable_index,
    dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the laplacian of the specified vector field.
   */
  vectorValue
  get_vector_laplacian(uint           global_variable_index,
                       dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the divergence of the specified vector field.
   */
  scalarValue
  get_vector_divergence(uint           global_variable_index,
                        dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the symmetric gradient of the specified vector field.
   */
  vectorGrad
  get_vector_symmetric_gradient(
    uint           global_variable_index,
    dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Return the curl of the specified vector field. Note that this is scalarValue
   * type for 2D and vectorValue type for 3D.
   */
  dealii::Tensor<1, (dim == 2 ? 1 : dim), dealii::VectorizedArray<number>>
  get_vector_curl(uint           global_variable_index,
                  dependencyType dependency_type = dependencyType::NORMAL) const;

  /**
   * \brief Set the residual value of the specified scalar field.
   */
  void
  set_scalar_value_term(uint           global_variable_index,
                        scalarValue    val,
                        dependencyType dependency_type = dependencyType::NORMAL);

  /**
   * \brief Set the residual gradient of the specified scalar field.
   */
  void
  set_scalar_gradient_term(uint           global_variable_index,
                           scalarGrad     grad,
                           dependencyType dependency_type = dependencyType::NORMAL);

  /**
   * \brief Set the residual value of the specified vector field.
   */
  void
  set_vector_value_term(uint           global_variable_index,
                        vectorValue    val,
                        dependencyType dependency_type = dependencyType::NORMAL);

  /**
   * \brief Set the residual gradient of the specified vector field.
   */
  void
  set_vector_gradient_term(uint           global_variable_index,
                           vectorGrad     grad,
                           dependencyType dependency_type = dependencyType::NORMAL);

  /**
   * \brief TODO: Add comments
   */
  void
  eval_local_operator(
    const std::function<void(variableContainer &,
                             const dealii::Point<dim, dealii::VectorizedArray<number>> &)>
                                                       &func,
    dealii::LinearAlgebra::distributed::Vector<number> &dst,
    const LinearAlgebra::distributed::Vector<number>   &src,
    const std::pair<uint, uint>                        &cell_range);

  /**
   * \brief TODO: Add comments
   */
  void
  eval_local_diagonal(
    const std::function<void(variableContainer &,
                             const dealii::Point<dim, dealii::VectorizedArray<number>> &)>
                                                       &func,
    dealii::LinearAlgebra::distributed::Vector<number> &dst,
    const std::pair<uint, uint>                        &cell_range,
    const uint                                         &global_var_index);

private:
  using scalar_FEEval = dealii::FEEvaluation<dim, degree, degree + 1, 1, number>;
  using vector_FEEval = dealii::FEEvaluation<dim, degree, degree + 1, dim, number>;

  /**
   * \brief Check that a variable value/gradient/hessians was marked as needed and thus
   * properly initialized.
   */
  void
  access_valid(const uint           &global_variable_index,
               const dependencyType &dependency_type,
               const EvalFlags      &flag) const;

  /**
   * \brief Return the number of quadrature points.
   */
  [[nodiscard]] uint
  get_n_q_points() const;

  /**
   * \brief Return the quadrate point location.
   */
  [[nodiscard]] dealii::Point<dim, dealii::VectorizedArray<number>>
  get_q_point_location() const;

  /**
   * \brief Initialize, read DOFs, and set evaulation flags for each variable.
   */
  void
  reinit_and_eval(const dealii::LinearAlgebra::distributed::Vector<number> &src,
                  uint                                                      cell);

  /**
   * \brief Initialize the cell for all dependencies of a certain variable index.
   */
  void
  reinit(uint cell, const uint &global_variable_index);

  /**
   * \brief Evaluate the flags on the cell for all dependencies of a certain variable
   * index.
   */
  void
  eval(const uint &global_variable_index);

  /**
   * \brief Integrate the residuals and distribute from local to global.
   */
  void
  integrate_and_distribute(dealii::LinearAlgebra::distributed::Vector<number> &dst);

  /**
   * \brief Integrate the residuals for a certain variable index.
   */
  void
  integrate(const uint &global_variable_index);

  /**
   * \brief Map of FEEvaluation objects for each active scalar variables. The first
   * mapping is for the global variable and the second is for the dependencyType.
   */
  std::unordered_map<uint,
                     std::unordered_map<dependencyType, std::unique_ptr<scalar_FEEval>>>
    scalar_vars_map;

  /**
   * \brief Map of FEEvaluation objects for each active vector variables. The first
   * mapping is for the global variable and the second is for the dependencyType.
   */
  std::unordered_map<uint,
                     std::unordered_map<dependencyType, std::unique_ptr<vector_FEEval>>>
    vector_vars_map;

  /**
   * \brief The attribute list of the relevant subset of variables.
   */
  const AttributesList &subset_attributes;

  /**
   * \brief The solve type
   */
  const solveType solve_type;

  /**
   * \brief The number of variables.
   */
  uint n_variables;

  /**
   * \brief The quadrature point index.
   */
  uint q_point = 0;

  /**
   * \brief Number of DoFs per cell.
   */
  uint n_dofs_per_cell;

  /**
   * \brief Diagonal matrix that is used for preconditioning.
   */
  std::unique_ptr<dealii::AlignedVector<dealii::VectorizedArray<number>>> diagonal;
};

#endif
