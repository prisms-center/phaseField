// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/vectorization.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <prismspf/core/field_container.h>
#include <prismspf/core/solution_indexer.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/group_solver_base.h>

#include <prismspf/config.h>

#include "prismspf/core/field_attributes.h"
#include "prismspf/core/pde_operator.h"
#include "prismspf/core/solve_group.h"

#if DEAL_II_VERSION_MAJOR >= 9 && DEAL_II_VERSION_MINOR >= 7
#  include <deal.II/base/enable_observer_pointer.h>
#  define MATRIX_FREE_OPERATOR_BASE dealii::EnableObserverPointer
#else
#  include <deal.II/base/subscriptor.h>
#  define MATRIX_FREE_OPERATOR_BASE dealii::Subscriptor
#endif

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class exists to evaluate a single user-defined operator for the matrix-free
 * implementation of some PDE.
 * @note Information such as the pde operator are passed in at construction/initialization
 * rather than being passed in during function calls because in certain contexts (eg.
 * GMG,) the operator gets called by a dealii function, rather than by prismspf, so it
 * needs to act as a standalone operator.
 *
 * @tparam dim The number of dimensions in the problem.
 * @tparam degree The polynomial degree of the shape functions.
 * @tparam number Datatype to use for `LinearAlgebra::distributed::Vector<number>`. Either
 * double or float.
 */
template <unsigned int dim, unsigned int degree, typename number>
class MFOperator : public MATRIX_FREE_OPERATOR_BASE
{
public:
  using SolutionVector = SolutionIndexer<dim, number>::SolutionVector;
  using ScalarValue    = dealii::VectorizedArray<number>;

  /**
   * @brief Constructor.
   */
  explicit MFOperator(SolveGroup _solve_group);

  /**
   * @brief Initialize operator.
   */
  void
  initialize(const dealii::MatrixFree<dim, number, ScalarValue> &_data,
             const ElementVolume<dim, degree, number>           &_element_volume_handler,
             const std::vector<unsigned int>                    &selected_field_indexes =
               std::vector<unsigned int>());

  // public:
  /**
   * @brief Calls cell_loop on function that calls user-defiend operator
   */
  void
  compute_operator(SolutionVector &dst, const SolutionVector &src) const;

private:
  /**
   * @brief Calls user-defiend operator
   */
  void
  compute_local_operator(const dealii::MatrixFree<dim, number, ScalarValue> &_data,
                         SolutionVector                                     &dst,
                         const SolutionVector                               &src,
                         const std::pair<unsigned int, unsigned int> &cell_range) const;

public:
  /**
   * @brief Compute the diagonal of this operator.
   */
  void
  compute_diagonal();

private:
  /**
   * @brief Local computation of the diagonal of the operator.
   */
  void
  compute_local_diagonal(const dealii::MatrixFree<dim, number, ScalarValue> &_data,
                         SolutionVector                                     &dst,
                         const unsigned int                                 &dummy,
                         const std::pair<unsigned int, unsigned int> &cell_range) const;

public:
  /**
   * @brief Compute the element volume. (And store in this object?)
   */
  void
  compute_element_volume();

private:
  /**
   * @brief Local computation of the element volume.
   */
  void
  compute_local_element_volume(
    const dealii::MatrixFree<dim, number, ScalarValue> &_data,
    SolutionVector                                     &dst,
    const unsigned int                                 &dummy,
    const std::pair<unsigned int, unsigned int>        &cell_range) const;

public:
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
  initialize_dof_vector(SolutionVector &dst, unsigned int dof_handler_index = 0) const;

  /**
   * @brief Set constrained entries to one.
   */
  void
  set_constrained_entries_to_one(SolutionVector &dst) const;

  /**
   * @brief Get read access to the MatrixFree object stored with this operator.
   */
  std::shared_ptr<const dealii::MatrixFree<dim, number, ScalarValue>>
  get_matrix_free() const;

  /**
   * @brief Get read access to the inverse diagonal of this operator.
   */
  const std::shared_ptr<dealii::DiagonalMatrix<SolutionVector>> &
  get_matrix_diagonal_inverse() const;

  /**
   * @brief Matrix-vector multiplication.
   */
  void
  vmult(SolutionVector &dst, const SolutionVector &src) const;

  // NOLINTBEGIN(readability-identifier-naming)

  /**
   * @brief Transpose matrix-vector multiplication.
   */
  void
  Tvmult(SolutionVector &dst, const SolutionVector &src) const;

  // NOLINTEND(readability-identifier-naming)

private:
  /**
   * @brief The attribute list of the relevant variables.
   */
  std::vector<FieldAttributes> field_attributes;

  /**
   * @brief The group being solved
   */
  SolveGroup solve_group;

  using Operator = void (PDEOperator<dim, degree, number>::*)(
    FieldContainer<dim, degree, number> &,                       /* variable_list */
    const dealii::Point<dim, dealii::VectorizedArray<number>> &, /* q_point_loc */
    const dealii::VectorizedArray<number> &,                     /* element_volume */
    unsigned int                                                 /* solve_group_id */
  ) const;
  /**
   * @brief PDE operator object (owning class instance of pde_op) for user defined PDEs.
   */
  const PDEOperator<dim, degree, number> *pde_operator;
  /**
   * @brief The actual PDE operator function ptr (eg. compute_rhs) for user defined PDEs.
   */
  Operator pde_op;

  /**
   * @brief Matrix-free object.
   */
  std::shared_ptr<const dealii::MatrixFree<dim, number, ScalarValue>> data;

  /**
   * @brief Indices of DoFs on edge in case the operator is used in GMG context.
   */
  std::vector<std::vector<unsigned int>> edge_constrained_indices;

  /**
   * @brief The diagonal matrix.
   */
  std::shared_ptr<dealii::DiagonalMatrix<SolutionVector>> diagonal_entries;

  /**
   * @brief The inverse diagonal matrix.
   */
  std::shared_ptr<dealii::DiagonalMatrix<SolutionVector>> inverse_diagonal_entries;
};

PRISMS_PF_END_NAMESPACE
