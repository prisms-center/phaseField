// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/vectorization.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <prismspf/core/dst_container.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/field_container.h>
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_indexer.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/group_solver_base.h>

#include <prismspf/config.h>

#include <vector>

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
  using BlockVector    = SolutionIndexer<dim, number>::BlockVector;
  using SolutionVector = SolutionIndexer<dim, number>::SolutionVector;
  using ScalarValue    = dealii::VectorizedArray<number>;
  using VectorValue    = dealii::Tensor<1, dim, ScalarValue>;
  using MatrixFree     = dealii::MatrixFree<dim, number, ScalarValue>;

  using Operator = void (PDEOperator<dim, degree, number>::*)(
    FieldContainer<dim, degree, number> &,                       /* variable_list */
    const dealii::Point<dim, dealii::VectorizedArray<number>> &, /* q_point_loc */
    const dealii::VectorizedArray<number> &,                     /* element_volume */
    unsigned int                                                 /* solve_group_id */
  ) const;

  using TensorRank = FieldInfo::TensorRank;
  template <TensorRank Rank>
  using Value = std::conditional_t<Rank == TensorRank::Scalar,
                                   ScalarValue,
                                   dealii::Tensor<int(Rank), dim, ScalarValue>>;

  template <TensorRank Rank>
  Value<Rank>
  identity()
  {
    static Value<Rank> ident = []()
    {
      Value<Rank> obj;
      for (int i = 0; i < Value<Rank>::n_independent_components; ++i)
        {
          obj[Value<Rank>::unrolled_to_component_indices(i)] = 1.0;
        }
    }();
    return ident;
  }

  template <>
  ScalarValue
  identity()
  {
    static ScalarValue ident(1.0);
    return ident;
  }

  template <TensorRank Rank>
  Value<Rank>
  zero()
  {
    return Value<Rank>();
  }

  template <>
  ScalarValue
  zero()
  {
    static ScalarValue zeroo(0.0);
    return zeroo;
  }

  /**
   * @brief Constructor.
   */
  explicit MFOperator(PDEOperator<dim, degree, number>   &operator_owner,
                      Operator                            oper,
                      const std::vector<FieldAttributes> &_field_attributes,
                      const SolutionIndexer<dim, number> &_solution_indexer,
                      unsigned int                        _relative_level,
                      DependencySet                       _dependency_map)
    : MATRIX_FREE_OPERATOR_BASE()
    , pde_operator(&operator_owner)
    , pde_op(oper)
    , field_attributes(_field_attributes)
    , solution_indexer(&_solution_indexer)
    , relative_level(_relative_level)
    , dependency_map(std::move(_dependency_map))
  {}

  /**
   * @brief Initialize.
   * @note This will look stylistically strange. We pass in the solution handler for the
   * fields being solved here, but we don't actually use any of the field data from it.
   * This is because the purpose of this MFOperator class is to provide an operator with
   * the signature `vmult(VectorType &dst, const VectorType &src)`. Because we are
   * using block vectors, for the MFOperator to work, it needs to have access to the index
   * mapping in addition to the matrix_free object. Both are stored in the solution
   * handler.
   */
  void
  initialize(const GroupSolutionHandler<dim, number> &dst_solution)
  {
    solve_group          = dst_solution.get_solve_group();
    data                 = dst_solution.get_matrix_free();
    field_to_block_index = dst_solution.get_global_to_block_index();
  }

  // public:
  /**
   * @brief Calls cell_loop on function that calls user-defined operator
   */
  void
  compute_operator(BlockVector &dst, const BlockVector &src) const;

private:
  /**
   * @brief Calls user-defiend operator
   */
  void
  compute_local_operator(const MatrixFree                            &_data,
                         BlockVector                                 &dst,
                         const BlockVector                           &src,
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
  compute_local_diagonal(const MatrixFree                            &_data,
                         BlockVector                                 &dst,
                         const unsigned int                          &dummy,
                         const std::pair<unsigned int, unsigned int> &cell_range) const;

  template <TensorRank Rank>
  dealii::AlignedVector<Value<Rank>>
  compute_field_diagonal(FieldContainer<dim, degree, number> &variable_list,
                         DSTContainer<dim, degree, number>   &dst_fields,
                         unsigned int                         field_index) const;

public:
  /**
   * @brief Compute the element volume. (And store in this object? void?)
   */
  void
  compute_element_volume();

private:
  /**
   * @brief Local computation of the element volume.
   */
  void
  compute_local_element_volume(
    const MatrixFree                            &_data,
    SolutionVector                              &dst,
    const unsigned int                          &dummy,
    const std::pair<unsigned int, unsigned int> &cell_range) const;

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
  std::shared_ptr<const MatrixFree>
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

  /**
   * @brief PDE operator object (owning class instance of pde_op) for user defined PDEs.
   */
  const PDEOperator<dim, degree, number> *pde_operator;
  /**
   * @brief The actual PDE operator function ptr (eg. compute_rhs) for user defined PDEs.
   */
  Operator pde_op;

  /**
   * @brief Read-access to fields.
   */
  const SolutionIndexer<dim, number> *solution_indexer;
  /**
   * @brief Level so that correct fields are read from indexer.
   */
  unsigned int relative_level;

  /**
   * @brief Which fields should be available to the solve.
   */
  DependencySet dependency_map;

  /**
   * @brief Mapping from field index to block index (only for dst).
   */
  std::vector<unsigned int> field_to_block_index;

  /**
   * @brief Matrix-free object.
   */
  std::shared_ptr<const MatrixFree> data;

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
