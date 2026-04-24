// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/field_container.h>
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/pde_operator_base.h>
#include <prismspf/core/simulation_timer.h>
#include <prismspf/core/solution_indexer.h>
#include <prismspf/core/solve_block.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/solver_base.h>

#include <prismspf/utilities/utilities.h>

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
  using ScalarValue = dealii::VectorizedArray<number>;
  using VectorValue = dealii::Tensor<1, dim, ScalarValue>;

  using Operator = void (PDEOperatorBase<dim, degree, number>::*)(
    FieldContainer<dim, degree, number> &, /* variable_list */
    const SimulationTimer &,               /* sim_timer */
    unsigned int                           /* solve_block_id */
  ) const;

  template <TensorRank Rank>
  using Value = std::conditional_t<Rank == TensorRank::Scalar,
                                   ScalarValue,
                                   dealii::Tensor<int(Rank), dim, ScalarValue>>;

  template <TensorRank Rank>
  static Value<Rank>
  identity()
  {
    if constexpr (Rank == TensorRank::Scalar)
      {
        return ScalarValue(1.0);
      }
    else
      {
        static Value<Rank> ident = []()
        {
          Value<Rank> obj;
          for (int i = 0; i < Value<Rank>::n_independent_components; ++i)
            {
              obj[Value<Rank>::unrolled_to_component_indices(i)] = 1.0;
            }
          return obj;
        }();
        return ident;
      }
  }

  template <TensorRank Rank>
  static Value<Rank>
  zero()
  {
    if constexpr (Rank == TensorRank::Scalar)
      {
        return ScalarValue(0.0);
      }
    else
      {
        return Value<Rank>();
      }
  }

  /**
   * @brief Constructor.
   * @note It might be better to provide a SolveContext object instead of individual
   * components
   */
  explicit MFOperator(const PDEOperatorBase<dim, degree, number> &operator_owner,
                      Operator                                    oper,
                      const std::vector<FieldAttributes>         &_field_attributes,
                      const SolutionIndexer<dim, number>         &_solution_indexer,
                      unsigned int                                _relative_level,
                      DependencyMap                               _dependency_map,
                      const SimulationTimer                      &_sim_timer)
    : MATRIX_FREE_OPERATOR_BASE()
    , pde_operator(&operator_owner)
    , pde_op(oper)
    , field_attributes(_field_attributes)
    , solution_indexer(&_solution_indexer)
    , relative_level(_relative_level)
    , dependency_map(std::move(_dependency_map))
    , sim_timer(&_sim_timer)
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
   *
   * Additionally, we modify dependency_map to include an entry for every field being
   * solved. This is to avoid a troublesome problem downstream in FieldContainer, where we
   * only want to initialize FEEvals for what is needed. (There could be a better way
   * of handling this, but this is the least complicated for now.)
   */
  void
  initialize(const GroupSolutionHandler<dim, number> &dst_solution)
  {
    solve_block          = dst_solution.get_solve_block();
    data                 = &(dst_solution.get_matrix_free(relative_level));
    field_to_block_index = dst_solution.get_global_to_block_index();
    for (unsigned int field_index : solve_block.field_indices)
      {
        dependency_map[field_index]; // creates entry if not already present
      }
  }

  // public:
  /**
   * @brief Calls cell_loop on function that calls user-defined operator
   * @note requires dst is not ghosted
   */
  void
  compute_operator(BlockVector<number>       &dst,
                   const BlockVector<number> &src = BlockVector<number>()) const;

private:
  /**
   * @brief Calls user-defined operator
   * @note requires dst is not ghosted
   */
  void
  compute_local_operator(const MatrixFree<dim, number>               &_data,
                         BlockVector<number>                         &dst,
                         const BlockVector<number>                   &src,
                         const std::pair<unsigned int, unsigned int> &cell_range) const;

public:
  /**
   * @brief Compute the diagonal of this operator.
   */
  void
  compute_diagonal(BlockVector<number> &dst, const BlockVector<number> &src) const;

private:
  /**
   * @brief Local computation of the diagonal of the operator.
   */
  void
  compute_local_diagonal(const MatrixFree<dim, number>               &_data,
                         BlockVector<number>                         &diagonal,
                         const BlockVector<number>                   &dummy_src,
                         const std::pair<unsigned int, unsigned int> &cell_range) const;

  template <TensorRank Rank>
  dealii::AlignedVector<Value<Rank>>
  compute_local_field_diagonal(FieldContainer<dim, degree, number> &variable_list,
                               BlockVector<number>                 &diagonal,
                               unsigned int                         field_index) const;

public:
  /**
   * @brief Set scaling diagonal
   */
  void
  set_scaling_diagonal(bool                                               scale,
                       const std::vector<const SolutionVector<number> *> &diagonal)
  {
    scaling_diagonal  = diagonal;
    scale_by_diagonal = scale;
  }

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
   * @brief Set constrained entries to one.
   */
  // void
  // set_constrained_entries_to_one(SolutionVector<number> &dst) const;

  /**
   * @brief Get read access to the MatrixFree<dim, number> object stored with this
   * operator.
   */
  const MatrixFree<dim, number> *
  get_matrix_free() const;

  /**
   * @brief Get read access to the inverse diagonal of this operator.
   */
  const std::shared_ptr<dealii::DiagonalMatrix<SolutionVector<number>>> &
  get_matrix_diagonal_inverse() const;

  /**
   * @brief Matrix-vector multiplication.
   * @note requires dst is not ghosted
   */
  void
  vmult(BlockVector<number> &dst, const BlockVector<number> &src) const;

  // NOLINTBEGIN(readability-identifier-naming)

  /**
   * @brief Transpose matrix-vector multiplication.
   */
  void
  Tvmult(BlockVector<number> &dst, const BlockVector<number> &src) const;

  // NOLINTEND(readability-identifier-naming)

  /**
   * @brief Whether to read plain dof values from src, otherwise applies homogeneous part
   * of constraints to the read of src.
   */
  bool read_plain = false;

private:
  /**
   * @brief The attribute list of the relevant variables.
   */
  std::vector<FieldAttributes> field_attributes;

  /**
   * @brief The block being solved
   */
  SolveBlock solve_block;

  /**
   * @brief PDE operator object (owning class instance of pde_op) for user defined PDEs.
   */
  const PDEOperatorBase<dim, degree, number> *pde_operator;
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
  DependencyMap dependency_map;

  /**
   * @brief Simulation timer
   */
  const SimulationTimer *sim_timer;

  /**
   * @brief Result of operator gets scaled by this (invm for explicit fields)
   */
  std::vector<const SolutionVector<number> *> scaling_diagonal;

  /**
   * @brief Whether or not to scale after operator result
   */
  bool scale_by_diagonal = false;

  /**
   * @brief Mapping from field index to block index (only for dst).
   */
  std::vector<unsigned int> field_to_block_index;

  /**
   * @brief Matrix-free object.
   */
  const MatrixFree<dim, number> *data;

  /**
   * @brief Indices of DoFs on edge in case the operator is used in GMG context.
   */
  std::vector<std::vector<unsigned int>> edge_constrained_indices;

  /**
   * @brief The diagonal matrix.
   */
  std::shared_ptr<dealii::DiagonalMatrix<SolutionVector<number>>> diagonal_entries;

  /**
   * @brief The inverse diagonal matrix.
   */
  std::shared_ptr<dealii::DiagonalMatrix<SolutionVector<number>>>
    inverse_diagonal_entries;
};

PRISMS_PF_END_NAMESPACE
