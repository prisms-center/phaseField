// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/vectorization.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#if DEAL_II_VERSION_MAJOR >= 9 && DEAL_II_VERSION_MINOR >= 7
#  include <deal.II/base/enable_observer_pointer.h>
#  define MATRIX_FREE_OPERATOR_BASE dealii::EnableObserverPointer
#else
#  include <deal.II/base/subscriptor.h>
#  define MATRIX_FREE_OPERATOR_BASE dealii::Subscriptor
#endif

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

struct VariableAttributes;

template <unsigned int dim, unsigned int degree, typename number>
class VariableContainer;

template <unsigned int dim, unsigned int degree, typename number>
class PDEOperator;

template <unsigned int dim, unsigned int degree, typename number>
class ElementVolume;

/**
 * @brief This is the abstract base class for the matrix-free implementation of some
 * PDE.
 *
 * @tparam dim The number of dimensions in the problem.
 * @tparam degree The polynomial degree of the shape functions.
 * @tparam number Datatype to use for `LinearAlgebra::distributed::Vector<number>`. Either
 * double or float.
 */
template <unsigned int dim, unsigned int degree, typename number>
class MatrixFreeOperator : public MATRIX_FREE_OPERATOR_BASE
{
public:
  using VectorType = dealii::LinearAlgebra::distributed::Vector<number>;
  using SizeType   = dealii::VectorizedArray<number>;

  /**
   * @brief Default constructor.
   *
   * TODO (landinjm): Should we have a default constructor and pass everything through
   * initialize? Need to pick one.
   */
  MatrixFreeOperator(
    std::map<unsigned int, VariableAttributes>              _attributes_list,
    std::shared_ptr<const PDEOperator<dim, degree, number>> _pde_operator,
    Types::Index                                            _solve_block,
    Types::Index _index             = Numbers::invalid_index,
    bool         _use_local_mapping = false);

  /**
   * @brief Initialize operator.
   */
  void
  initialize(std::shared_ptr<const dealii::MatrixFree<dim, number, SizeType>> _data,
             const ElementVolume<dim, degree, number> &_element_volume_handler,
             const std::vector<unsigned int>          &selected_field_indexes =
               std::vector<unsigned int>());

  /**
   * @brief Return the number of DoFs.
   */
  dealii::types::global_dof_index
  m() const;

  /**
   * @brief Return the value of the matrix entry. This function is only valid when row ==
   * col and when the diagonal is initialized. Additionally, this is only used so that we
   * may compile. Trying to use this function will throw an error.
   */
  number
  el(const unsigned int &row, const unsigned int &col) const;

  /**
   * @brief Release all memory and return to state like having called the default
   * constructor.
   */
  void
  clear();

  /**
   * @brief Initialize a given vector with the MatrixFree object that this object
   * contains.
   */
  void
  initialize_dof_vector(VectorType &dst, unsigned int dof_handler_index = 0) const;

  /**
   * @brief Set constrained entries to one.
   */
  void
  set_constrained_entries_to_one(VectorType &dst) const;

  /**
   * @brief Get read access to the MatrixFree object stored with this operator.
   */
  std::shared_ptr<const dealii::MatrixFree<dim, number, SizeType>>
  get_matrix_free() const;

  /**
   * @brief Get read access to the inverse diagonal of this operator.
   */
  const std::shared_ptr<dealii::DiagonalMatrix<VectorType>> &
  get_matrix_diagonal_inverse() const;

  // cppcheck-suppress-begin passedByValue

  /**
   * @brief Add the mappings from global to local solution vectors.
   */
  void
  add_global_to_local_mapping(std::vector<Types::Index> _global_to_local_solution);

  /**
   * @brief Add the solution subset for src vector.
   */
  void
  add_src_solution_subset(
    std::vector<VectorType *> _src_solution_subset = std::vector<VectorType *>());

  // cppcheck-suppress-end passedByValue

  /**
   * @brief Matrix-vector multiplication.
   */
  void
  vmult(VectorType &dst, const VectorType &src) const;

  // NOLINTBEGIN(readability-identifier-naming)

  /**
   * @brief Transpose matrix-vector multiplication.
   */
  void
  Tvmult(VectorType &dst, const VectorType &src) const;

  // NOLINTEND(readability-identifier-naming)

  /**
   * @brief Compute the explicit update.
   */
  void
  compute_explicit_update(std::vector<VectorType *>       &dst,
                          const std::vector<VectorType *> &src) const;

  /**
   * @brief Compute the explicit update for postprocessed fields.
   */
  void
  compute_postprocess_explicit_update(std::vector<VectorType *>       &dst,
                                      const std::vector<VectorType *> &src) const;

  /**
   * @brief Compute a nonexplicit auxiliary update.
   */
  void
  compute_nonexplicit_auxiliary_update(std::vector<VectorType *>       &dst,
                                       const std::vector<VectorType *> &src) const;

  /**
   * @brief Compute the residual of this operator. This is the b in Ax=b.
   */
  void
  compute_residual(VectorType &dst, const VectorType &src) const;

  /**
   * @brief Compute the diagonal of this operator.
   */
  void
  compute_diagonal(unsigned int field_index);

private:
  /**
   * @brief Local computation of the explicit update.
   */
  void
  compute_local_explicit_update(
    const dealii::MatrixFree<dim, number, SizeType> &_data,
    std::vector<VectorType *>                       &dst,
    const std::vector<VectorType *>                 &src,
    const std::pair<unsigned int, unsigned int>     &cell_range) const;

  /**
   * @brief Local computation of the explicit update of postprocessed fields.
   */
  void
  compute_local_postprocess_explicit_update(
    const dealii::MatrixFree<dim, number, SizeType> &_data,
    std::vector<VectorType *>                       &dst,
    const std::vector<VectorType *>                 &src,
    const std::pair<unsigned int, unsigned int>     &cell_range) const;

  /**
   * @brief Local computation of the nonexplicit auxiliary update.
   */
  void
  compute_local_nonexplicit_auxiliary_update(
    const dealii::MatrixFree<dim, number, SizeType> &_data,
    std::vector<VectorType *>                       &dst,
    const std::vector<VectorType *>                 &src,
    const std::pair<unsigned int, unsigned int>     &cell_range) const;

  /**
   * @brief Local computation of the residual of the operator.
   */
  void
  compute_local_residual(const dealii::MatrixFree<dim, number, SizeType> &_data,
                         VectorType                                      &dst,
                         const VectorType                                &src,
                         const std::pair<unsigned int, unsigned int> &cell_range) const;

  /**
   * @brief Local computation of the newton update of the operator.
   */
  void
  compute_local_newton_update(
    const dealii::MatrixFree<dim, number, SizeType> &_data,
    VectorType                                      &dst,
    const VectorType                                &src,
    const std::pair<unsigned int, unsigned int>     &cell_range) const;

  /**
   * @brief Local computation of the diagonal of the operator.
   */
  void
  local_compute_diagonal(const dealii::MatrixFree<dim, number, SizeType> &_data,
                         VectorType                                      &dst,
                         const unsigned int                              &dummy,
                         const std::pair<unsigned int, unsigned int> &cell_range) const;

  /**
   * @brief The attribute list of the relevant variables.
   */
  std::map<unsigned int, VariableAttributes> attributes_list;

  /**
   * @brief The element volume container
   */
  const ElementVolume<dim, degree, number> *element_volume_handler = nullptr;

  /**
   * @brief PDE operator object for user defined PDEs.
   */
  std::shared_ptr<const PDEOperator<dim, degree, number>> pde_operator;

  /**
   * @brief The solve block that is being evaluated
   */
  Types::Index solve_block = 0;

  /**
   * @brief Field index that is being evaluated.
   */
  Types::Index index = Numbers::invalid_index;

  /**
   * @brief Whether to use local mapping for the VariableContainer object.
   */
  bool use_local_mapping = false;

  /**
   * @brief Matrix-free object.
   */
  std::shared_ptr<const dealii::MatrixFree<dim, number, SizeType>> data;

  /**
   * @brief Selected fields for which we'll evaluate.
   */
  std::vector<unsigned int> selected_fields;

  /**
   * @brief Indices of DoFs on edge in case the operator is used in GMG context.
   */
  std::vector<std::vector<unsigned int>> edge_constrained_indices;

  /**
   * @brief Mapping from global solution vectors to the local ones
   */
  std::vector<Types::Index> global_to_local_solution;

  /**
   * @brief Subset of fields that are necessary for the source.
   */
  std::vector<VectorType *> src_solution_subset;

  /**
   * @brief The diagonal matrix.
   */
  std::shared_ptr<dealii::DiagonalMatrix<VectorType>> diagonal_entries;

  /**
   * @brief The inverse diagonal matrix.
   */
  std::shared_ptr<dealii::DiagonalMatrix<VectorType>> inverse_diagonal_entries;
};

PRISMS_PF_END_NAMESPACE
