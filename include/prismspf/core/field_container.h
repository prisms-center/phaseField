// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/lac/vector.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/dependencies.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solution_indexer.h>
#include <prismspf/core/solve_block.h>
#include <prismspf/core/system_wide.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class permits the access of a subset of indexed fields and gives an error
 * if any non-allowed fields are requested.
 *
 * @tparam dim The number of dimensions in the problem.
 * @tparam degree The polynomial degree of the shape functions.
 * @tparam number Datatype to use for `dealii::VectorizedArray<number>`. Either
 * double or float.
 *
 * Importantly, this class is mostly a wrapper for the dealii::FEEvaluation class, which
 * allows for the evaluation of functions at quadrature points and cell integrations
 * (basically the backbone of the FEM solver). As such, it is one of the most important
 * classes when it comes to performance because these functions are called millions of
 * times.
 */
template <unsigned int dim, unsigned int degree, typename number>
class FieldContainer
{
public:
  /**
   * @brief Typedef for the basic value that the user manipulates.
   */
  using ScalarValue = dealii::VectorizedArray<number>;

  /**
   * @brief Typedef for the basic vector that the user manipulates.
   */
  using VectorValue = dealii::Tensor<1, dim, ScalarValue>;

  template <TensorRank Rank>
  using Value = std::conditional_t<Rank == TensorRank::Scalar,
                                   ScalarValue,
                                   dealii::Tensor<int(Rank), dim, ScalarValue>>;

  template <TensorRank Rank>
  using Gradient = dealii::Tensor<int(Rank) + 1, dim, ScalarValue>;

  template <TensorRank Rank>
  using Hessian = dealii::Tensor<int(Rank) + 2, dim, ScalarValue>;

  template <TensorRank Rank>
  using FEEval =
    dealii::FEEvaluation<dim,
                         degree,
                         degree + 1,
                         dealii::Tensor<int(Rank), dim>::n_independent_components,
                         number,
                         ScalarValue>;

  /**
   * @brief Return the tensor rank from the specified template value.
   */
  template <typename ValType>
  static constexpr TensorRank RankFromVal = []() constexpr
  {
    if constexpr (std::is_same_v<ValType, ScalarValue> ||
                  std::is_same_v<ValType, typename ScalarValue::value_type>)
      {
        return TensorRank::Scalar;
      }
    else
      {
        return TensorRank(ValType::rank);
      }
  }();

  /**
   * @brief Return the tensor rank from the specified template gradient.
   */
  template <typename GradType>
  static constexpr TensorRank RankFromGrad = []() constexpr
  {
    if constexpr (std::is_same_v<GradType, ScalarValue>)
      {
        return TensorRank::Scalar;
      }
    else
      {
        return TensorRank(GradType::rank - 1);
      }
  }();

  /**
   * @brief Struct to hold the relevant dealii::FEEvaluation for a given solution block
   * index.
   */
  template <TensorRank Rank>
  struct FEEValuationDeps
  {
    /**
     * @brief dealii::FEEvaluation object and the evaluation flags.
     */
    using FEEDepPair = std::pair<FEEval<Rank>, EvalFlags>;

    /**
     * @brief Pointer the the FEEDepPair. We have to use ptrs because we don't want to
     * construct fields we don't need.
     */
    using FEEDepPairPtr = std::shared_ptr<FEEDepPair>;

    /**
     * @brief FEEvaluation for the current solution.
     */
    FEEDepPairPtr fe_eval;

    /**
     * @brief FEEvaluation for the src and dst solutions.
     */
    FEEDepPairPtr fe_eval_src_dst;

    /**
     * @brief Collection of FEEvaluation for old solutions (< n -1 state).
     */
    std::vector<FEEDepPairPtr> fe_eval_old;

    /**
     * @brief Destination integration flags. We use this to integrate value and gradient
     * terms, if they exist.
     */
    EvalFlags integration_flags = EvalFlags::nothing;

    /**
     * @brief The solution block.
     *
     * @remark It would look nicer to just use the SolutionIndexer, but this way decreases
     * indexing, increasing performance.
     */
    const SolutionLevel<dim, number> *solution_level = nullptr;

    /**
     * @brief The solution block index.
     *
     * This is the index that tells dealii::FEEvaluation the corresponding
     * dealii::DoFHandler, dealii::AffineConstraints, and dealii::Quadrature to use from
     * the dealii::MatrixFree.
     */
    unsigned int block_index = -1;

    /**
     * @brief Default constructor.
     */
    FEEValuationDeps() = default;

    FEEValuationDeps(
      const Dependency                                                  &dependency,
      const std::pair<const SolutionLevel<dim, number> *, unsigned int> &mf_id_pair,
      bool                                                               is_dst);

    template <DependencyType type>
    const FEEval<Rank> &
    get() const;

    template <DependencyType type>
    FEEval<Rank> &
    get();

    const FEEval<Rank> &
    get(DependencyType type) const;

    FEEval<Rank> &
    get(DependencyType type);

    void
    reinit(unsigned int cell);

    void
    eval(const BlockVector<number> *_src_solutions, bool plain);

    void
    reinit_and_eval(unsigned int               cell,
                    const BlockVector<number> *_src_solutions,
                    bool                       plain);

    void
    integrate();

    void
    distribute(BlockVector<number> *dst_solutions);

    void
    integrate_and_distribute(BlockVector<number> *dst_solutions);
  };

  /**
   * @brief Constructor.
   *
   * `field_attributes` is the collection of field attributes. We only use this for the
   * rank of field (Scalar/Vector).
   * `solution_indexer` is the solution indexer object, which allows us to get access to
   * solution vectors from in and outside the current solve block.
   * `relative_level` is the multigrid level.
   * `dependency_map` is the map of field indices and their dependency flags.
   * `solve_block` is the object that contains attributes about the current solve block.
   * `matrix_free` is the dealii::MatrixFree object.
   */
  FieldContainer(const std::vector<FieldAttributes> &_field_attributes,
                 const SolutionIndexer<dim, number> &_solution_indexer,
                 unsigned int                        _relative_level,
                 const DependencyMap                &dependency_map,
                 const SolveBlock                   &_solve_block,
                 const MatrixFree<dim, number>      &matrix_free);

  /**
   * @brief Initialize based on cell for all dependencies.
   */
  void
  reinit(unsigned int cell);

  /**
   * @brief Read solution vector, and evaluate based on dependency flags for all
   * dependencies.
   */
  void
  eval(const BlockVector<number> *src_solutions, bool plain);

  /**
   * @brief Initialize based on cell, read solution vector, and evaluate based on
   * dependency flags for all dependencies.
   *
   * This is more efficient than calling `reinit` and `eval` individually.
   */
  void
  reinit_and_eval(unsigned int               cell,
                  const BlockVector<number> *src_solutions,
                  bool                       plain);

  /**
   * @brief Integrate the residuals.
   */
  void
  integrate();

  /**
   * @brief Distribute the integrated residuals.
   */
  void
  distribute(BlockVector<number> *dst_solutions);

  /**
   * @brief Integrate the residuals and distribute from local to global.
   *
   * This is more efficient that calling `integrate` and `distribute` individually.
   */
  void
  integrate_and_distribute(BlockVector<number> *dst_solutions);

  /**
   * @brief Set the current quadrature point.
   */
  void
  set_q_point(unsigned int q);

  /**
   * @brief Return the value of the specified field.
   */
  template <TensorRank Rank, DependencyType type>
  [[nodiscard]] Value<Rank>
  get_value(Types::Index global_variable_index) const;

  /**
   * @brief Return the value of the specified field.
   */
  template <TensorRank Rank>
  [[nodiscard]] Value<Rank>
  get_value(Types::Index global_variable_index, DependencyType type) const;

  /**
   * @brief Return the gradient of the specified field.
   */
  template <TensorRank Rank, DependencyType type>
  [[nodiscard]] Gradient<Rank>
  get_gradient(Types::Index global_variable_index) const;

  /**
   * @brief Return the gradient of the specified field.
   */
  template <TensorRank Rank>
  [[nodiscard]] Gradient<Rank>
  get_gradient(Types::Index global_variable_index, DependencyType type) const;

  /**
   * @brief Return the hessian of the specified field.
   */
  template <TensorRank Rank, DependencyType type>
  [[nodiscard]] Hessian<Rank>
  get_hessian(Types::Index global_variable_index) const;

  /**
   * @brief Return the hessian of the specified field.
   */
  template <TensorRank Rank>
  [[nodiscard]] Hessian<Rank>
  get_hessian(Types::Index global_variable_index, DependencyType type) const;

  /**
   * @brief Return the diagonal of the hessian of the specified field.
   */
  template <TensorRank Rank, DependencyType type>
  [[nodiscard]] Gradient<Rank>
  get_hessian_diagonal(Types::Index global_variable_index) const;

  /**
   * @brief Return the diagonal of the hessian of the specified field.
   */
  template <TensorRank Rank>
  [[nodiscard]] Gradient<Rank>
  get_hessian_diagonal(Types::Index global_variable_index, DependencyType type) const;

  /**
   * @brief Return the laplacian of the specified field.
   */
  template <TensorRank Rank, DependencyType type>
  [[nodiscard]] Value<Rank>
  get_laplacian(Types::Index global_variable_index) const;

  /**
   * @brief Return the laplacian of the specified field.
   */
  template <TensorRank Rank>
  [[nodiscard]] Value<Rank>
  get_laplacian(Types::Index global_variable_index, DependencyType type) const;

  /**
   * @brief Return the divergence of the specified field.
   */
  template <TensorRank Rank, DependencyType type>
  [[nodiscard]] ScalarValue
  get_divergence(Types::Index global_variable_index) const;

  /**
   * @brief Return the divergence of the specified field.
   */
  template <TensorRank Rank>
  [[nodiscard]] ScalarValue
  get_divergence(Types::Index global_variable_index, DependencyType type) const;

  /**
   * @brief Return the symmetric gradient of the specified field.
   */
  template <TensorRank Rank, DependencyType type>
  [[nodiscard]] dealii::SymmetricTensor<2, dim, ScalarValue>
  get_symmetric_gradient(Types::Index global_variable_index) const;

  /**
   * @brief Return the symmetric gradient of the specified field.
   */
  template <TensorRank Rank>
  [[nodiscard]] dealii::SymmetricTensor<2, dim, ScalarValue>
  get_symmetric_gradient(Types::Index global_variable_index, DependencyType type) const;

  /**
   * @brief Return the curl of the specified field.
   */
  template <TensorRank Rank, DependencyType type>
  [[nodiscard]] dealii::Tensor<1, (dim == 2 ? 1 : dim), ScalarValue>
  get_curl(Types::Index global_variable_index) const;

  /**
   * @brief Return the curl of the specified field.
   */
  template <TensorRank Rank>
  [[nodiscard]] dealii::Tensor<1, (dim == 2 ? 1 : dim), ScalarValue>
  get_curl(Types::Index global_variable_index, DependencyType type) const;

  /**
   * @brief Return the quadrature point location.
   */
  [[nodiscard]] dealii::Point<dim, ScalarValue>
  get_q_point_location() const;

  /**
   * @brief Return the quadrature point location.
   */
  [[nodiscard]] ScalarValue
  get_element_volume() const;

  /**
   * @brief Return the number of quadrature points.
   */
  [[nodiscard]] unsigned int
  get_n_q_points() const;

  /**
   * @brief Set the residual value of the specified scalar/vector field.
   */
  template <typename ValType>
  void
  set_value_term(Types::Index global_variable_index, const ValType &val);

  /**
   * @brief Set the residual gradient of the specified scalar/vector field.
   */
  template <typename GradType>
  void
  set_gradient_term(Types::Index global_variable_index, const GradType &val);

private:
  template <TensorRank Rank>
  std::vector<FEEValuationDeps<Rank>> &
  get_relevant_feeval_vector();

  template <TensorRank Rank>
  const std::vector<FEEValuationDeps<Rank>> &
  get_relevant_feeval_vector() const;

  /**
   * @brief Check whether the entry for the FEEvaluation is within the bounds of the
   * vector.
   */
  void
  feevaluation_size_valid(Types::Index field_index) const;

  /**
   * @brief Check whether the entry for the FEEvaluation is within the bounds of the
   * vector and not a nullptr.
   */
  void
  feevaluation_exists(Types::Index field_index, Types::Index dependency_index) const;

  /**
   * @brief Check that a variable value/gradient/hessians was marked as needed and thus
   * properly initialized.
   */
  void
  access_valid(Types::Index                             field_index,
               DependencyType                           dependency_type,
               dealii::EvaluationFlags::EvaluationFlags flag) const;

  /**
   * @brief Check that a value is valid for submission.
   */
  void
  submission_valid(Types::Index field_index, DependencyType dependency_type) const;

  /**
   * @brief Pointer to the vector of field attributes.
   */
  const std::vector<FieldAttributes> *field_attributes_ptr;

  /**
   * @brief Pointer to the solution indexer.
   */
  const SolutionIndexer<dim, number> *solution_indexer;

  /**
   * @brief Collection of FEEvaluation dependencies for scalar fields.
   */
  std::vector<FEEValuationDeps<TensorRank::Scalar>> feeval_deps_scalar;

  /**
   * @brief Collection of FEEvaluation dependencies for vector fields.
   */
  std::vector<FEEValuationDeps<TensorRank::Vector>> feeval_deps_vector;

  /**
   * @brief Solve block information.
   */
  const SolveBlock *solve_block;

  /**
   * @brief FEEvaluation object for generic cell operations.
   *
   * This is what we use when we call `get_q_point_location`, `get_element_volume`, and
   * `get_n_q_points()`.
   *
   * @note This is constructed with the 0th index of the MatrixFree object. If that is
   * not a scalar you'll run into issues.
   */
  FEEval<TensorRank::Scalar> shared_feeval_scalar;

  /**
   * @brief Multigrid level.
   */
  unsigned int relative_level;

  /**
   * @brief The quadrature point.
   */
  unsigned int q_point = 0;
};

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
FieldContainer<dim, degree, number>::FEEValuationDeps<Rank>::FEEValuationDeps(
  const Dependency                                                  &dependency,
  const std::pair<const SolutionLevel<dim, number> *, unsigned int> &mf_id_pair,
  bool                                                               is_dst)
  : solution_level(mf_id_pair.first)
  , block_index(mf_id_pair.second)
{
  // Make an FEEvaluation if the current solution needs to be evaluated
  if (dependency.flag)
    {
      fe_eval = std::make_shared<FEEDepPair>(FEEval<Rank>(solution_level->matrix_free,
                                                          block_index),
                                             dependency.flag);
    }
  // Make FEEvaluations for the the old solutions
  fe_eval_old.resize(dependency.old_flags.size(), nullptr);
  for (unsigned int age = 0; age < dependency.old_flags.size(); ++age)
    {
      if (dependency.old_flags.at(age)) // kinda redundant... maybe remove?
        {
          fe_eval_old[age] =
            std::make_shared<FEEDepPair>(FEEval<Rank>(solution_level->matrix_free,
                                                      block_index),
                                         dependency.old_flags.at(age));
        }
    }
  if (dependency.src_flag || is_dst)
    {
      fe_eval_src_dst =
        std::make_shared<FEEDepPair>(FEEval<Rank>(solution_level->matrix_free,
                                                  block_index),
                                     dependency.src_flag);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
template <DependencyType Type>
inline DEAL_II_ALWAYS_INLINE const typename FieldContainer<dim, degree, number>::
  template FEEval<Rank> &
  FieldContainer<dim, degree, number>::FEEValuationDeps<Rank>::get() const
{
  // TODO: Assertions
  if constexpr (Type == DependencyType::SRC || Type == DependencyType::DST)
    {
      return fe_eval_src_dst->first;
    }
  else if constexpr (Type == DependencyType::Current)
    {
      return fe_eval->first;
    }
  else
    {
      return fe_eval_old[int(Type) - 1]->first;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
template <DependencyType Type>
inline DEAL_II_ALWAYS_INLINE
  typename FieldContainer<dim, degree, number>::template FEEval<Rank> &
  FieldContainer<dim, degree, number>::FEEValuationDeps<Rank>::get()
{
  // TODO: Assertions
  if constexpr (Type == DependencyType::SRC || Type == DependencyType::DST)
    {
      return fe_eval_src_dst->first;
    }
  else if constexpr (Type == DependencyType::Current)
    {
      return fe_eval->first;
    }
  else
    {
      return fe_eval_old[int(Type) - 1]->first;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline DEAL_II_ALWAYS_INLINE const typename FieldContainer<dim, degree, number>::
  template FEEval<Rank> &
  FieldContainer<dim, degree, number>::FEEValuationDeps<Rank>::get(
    DependencyType type) const
{
  // TODO: Assertions
  if (type == DependencyType::SRC || type == DependencyType::DST)
    {
      return fe_eval_src_dst->first;
    }
  if (type == DependencyType::Current)
    {
      return fe_eval->first;
    }
  {
    return fe_eval_old[type - 1]->first;
  }
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline DEAL_II_ALWAYS_INLINE
  typename FieldContainer<dim, degree, number>::template FEEval<Rank> &
  FieldContainer<dim, degree, number>::FEEValuationDeps<Rank>::get(DependencyType type)
{
  // TODO: Assertions
  if (type == DependencyType::SRC || type == DependencyType::DST)
    {
      return fe_eval_src_dst->first;
    }
  if (type == DependencyType::Current)
    {
      return fe_eval->first;
    }
  {
    return fe_eval_old[type - 1]->first;
  }
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline void
FieldContainer<dim, degree, number>::FEEValuationDeps<Rank>::reinit(unsigned int cell)
{
  if (fe_eval)
    {
      fe_eval->first.reinit(cell);
    }
  for (auto &old_fe_eval : fe_eval_old)
    {
      if (old_fe_eval)
        {
          old_fe_eval->first.reinit(cell);
        }
    }
  if (fe_eval_src_dst)
    {
      fe_eval_src_dst->first.reinit(cell);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline void
FieldContainer<dim, degree, number>::FEEValuationDeps<Rank>::eval(
  const BlockVector<number> *_src_solutions,
  bool                       plain)
{
  // NOTE: `read_dof_values_plain` must be called here so that constraints aren't
  // implicitly applied. This allows us to have inhomogeneous constraints.
  if (fe_eval)
    {
      fe_eval->first.read_dof_values_plain(solution_level->solutions.block(block_index));
      fe_eval->first.evaluate(fe_eval->second);
    }
  for (unsigned int age = 0; age < fe_eval_old.size(); ++age)
    {
      if (FEEDepPairPtr &old_fe_eval = fe_eval_old[age])
        {
          old_fe_eval->first.read_dof_values_plain(
            solution_level->old_solutions[age].block(block_index));
          old_fe_eval->first.evaluate(old_fe_eval->second);
        }
    }
  if (fe_eval_src_dst && fe_eval_src_dst->second != EvalFlags::nothing)
    {
      if (plain)
        {
          fe_eval_src_dst->first.read_dof_values_plain(
            _src_solutions->block(block_index));
        }
      else
        {
          fe_eval_src_dst->first.read_dof_values(_src_solutions->block(block_index));
        }
      fe_eval_src_dst->first.evaluate(fe_eval_src_dst->second);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline void
FieldContainer<dim, degree, number>::FEEValuationDeps<Rank>::reinit_and_eval(
  unsigned int               cell,
  const BlockVector<number> *_src_solutions,
  bool                       plain)
{
  // NOTE: `read_dof_values_plain` must be called here so that constraints aren't
  // implicitly applied. This allows us to have inhomogeneous constraints.
  if (fe_eval)
    {
      fe_eval->first.reinit(cell);
      fe_eval->first.read_dof_values_plain(solution_level->solutions.block(block_index));
      fe_eval->first.evaluate(fe_eval->second);
    }
  for (unsigned int age = 0; age < fe_eval_old.size(); ++age)
    {
      if (FEEDepPairPtr &old_fe_eval = fe_eval_old[age])
        {
          old_fe_eval->first.reinit(cell);
          old_fe_eval->first.read_dof_values_plain(
            solution_level->old_solutions[age].block(block_index));
          old_fe_eval->first.evaluate(old_fe_eval->second);
        }
    }
  if (fe_eval_src_dst)
    {
      fe_eval_src_dst->first.reinit(cell);
      if (fe_eval_src_dst->second != EvalFlags::nothing)
        {
          if (plain)
            {
              fe_eval_src_dst->first.read_dof_values_plain(
                _src_solutions->block(block_index));
            }
          else
            {
              fe_eval_src_dst->first.read_dof_values(_src_solutions->block(block_index));
            }
          fe_eval_src_dst->first.evaluate(fe_eval_src_dst->second);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline void
FieldContainer<dim, degree, number>::FEEValuationDeps<Rank>::integrate()
{
  if (fe_eval_src_dst)
    {
      fe_eval_src_dst->first.integrate(integration_flags);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline void
FieldContainer<dim, degree, number>::FEEValuationDeps<Rank>::distribute(
  BlockVector<number> *dst_solutions)
{
  if (fe_eval_src_dst)
    {
      fe_eval_src_dst->first.distribute_local_to_global(
        dst_solutions->block(block_index));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline void
FieldContainer<dim, degree, number>::FEEValuationDeps<Rank>::integrate_and_distribute(
  BlockVector<number> *dst_solutions)
{
  if (fe_eval_src_dst)
    {
      fe_eval_src_dst->first.integrate_scatter(integration_flags,
                                               dst_solutions->block(block_index));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
inline void
FieldContainer<dim, degree, number>::reinit(unsigned int cell)
{
  for (auto &fe_eval : feeval_deps_scalar)
    {
      fe_eval.reinit(cell);
    }
  for (auto &fe_eval : feeval_deps_vector)
    {
      fe_eval.reinit(cell);
    }
  shared_feeval_scalar.reinit(cell);
}

template <unsigned int dim, unsigned int degree, typename number>
inline void
FieldContainer<dim, degree, number>::eval(const BlockVector<number> *src_solutions,
                                          bool                       plain)
{
  for (auto &fe_eval : feeval_deps_scalar)
    {
      fe_eval.eval(src_solutions, plain);
    }
  for (auto &fe_eval : feeval_deps_vector)
    {
      fe_eval.eval(src_solutions, plain);
    }
  // Don't eval `shared_feeval_scalar` because we only use it for information.
}

template <unsigned int dim, unsigned int degree, typename number>
inline void
FieldContainer<dim, degree, number>::reinit_and_eval(
  unsigned int               cell,
  const BlockVector<number> *src_solutions,
  bool                       plain)
{
  for (auto &fe_eval : feeval_deps_scalar)
    {
      fe_eval.reinit_and_eval(cell, src_solutions, plain);
    }
  for (auto &fe_eval : feeval_deps_vector)
    {
      fe_eval.reinit_and_eval(cell, src_solutions, plain);
    }
  // Don't eval `shared_feeval_scalar` because we only use it for information.
  shared_feeval_scalar.reinit(cell);
}

template <unsigned int dim, unsigned int degree, typename number>
inline void
FieldContainer<dim, degree, number>::integrate()
{
  const std::vector<FieldAttributes> &field_attributes = *field_attributes_ptr;
  for (const Types::Index &field_index : solve_block->field_indices)
    {
      if (field_attributes[field_index].field_type == TensorRank::Scalar)
        {
          feeval_deps_scalar[field_index].integrate();
        }
      else /* vector */
        {
          feeval_deps_vector[field_index].integrate();
        }
    }
  // Don't integrate `shared_feeval_scalar` because we only use it for information.
}

template <unsigned int dim, unsigned int degree, typename number>
inline void
FieldContainer<dim, degree, number>::distribute(BlockVector<number> *dst_solutions)
{
  const std::vector<FieldAttributes> &field_attributes = *field_attributes_ptr;
  for (const Types::Index &field_index : solve_block->field_indices)
    {
      if (field_attributes[field_index].field_type == TensorRank::Scalar)
        {
          feeval_deps_scalar[field_index].distribute(dst_solutions);
        }
      else /* vector */
        {
          feeval_deps_vector[field_index].distribute(dst_solutions);
        }
    }
  // Don't distribute `shared_feeval_scalar` because we only use it for information.
}

template <unsigned int dim, unsigned int degree, typename number>
inline void
FieldContainer<dim, degree, number>::integrate_and_distribute(
  BlockVector<number> *dst_solutions)
{
  const std::vector<FieldAttributes> &field_attributes = *field_attributes_ptr;
  for (const Types::Index &field_index : solve_block->field_indices)
    {
      if (field_attributes[field_index].field_type == TensorRank::Scalar)
        {
          feeval_deps_scalar[field_index].integrate_and_distribute(dst_solutions);
        }
      else /* vector */
        {
          feeval_deps_vector[field_index].integrate_and_distribute(dst_solutions);
        }
    }
  // Don't integrate and distribute `shared_feeval_scalar` because we only use it for
  // information.
}

template <unsigned int dim, unsigned int degree, typename number>
inline DEAL_II_ALWAYS_INLINE void
FieldContainer<dim, degree, number>::set_q_point(unsigned int q)
{
  q_point = q;
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank, DependencyType type>
inline DEAL_II_ALWAYS_INLINE
  typename FieldContainer<dim, degree, number>::template Value<Rank>
  FieldContainer<dim, degree, number>::get_value(Types::Index global_variable_index) const
{
  return get_relevant_feeval_vector<Rank>()[global_variable_index]
    .template get<type>()
    .get_value(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline DEAL_II_ALWAYS_INLINE
  typename FieldContainer<dim, degree, number>::template Value<Rank>
  FieldContainer<dim, degree, number>::get_value(Types::Index   global_variable_index,
                                                 DependencyType type) const
{
  return get_relevant_feeval_vector<Rank>()[global_variable_index].get(type).get_value(
    q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank, DependencyType type>
inline DEAL_II_ALWAYS_INLINE
  typename FieldContainer<dim, degree, number>::template Gradient<Rank>
  FieldContainer<dim, degree, number>::get_gradient(
    Types::Index global_variable_index) const
{
  return get_relevant_feeval_vector<Rank>()[global_variable_index]
    .template get<type>()
    .get_gradient(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline DEAL_II_ALWAYS_INLINE
  typename FieldContainer<dim, degree, number>::template Gradient<Rank>
  FieldContainer<dim, degree, number>::get_gradient(Types::Index   global_variable_index,
                                                    DependencyType type) const
{
  return get_relevant_feeval_vector<Rank>()[global_variable_index].get(type).get_gradient(
    q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank, DependencyType type>
inline DEAL_II_ALWAYS_INLINE
  typename FieldContainer<dim, degree, number>::template Hessian<Rank>
  FieldContainer<dim, degree, number>::get_hessian(
    Types::Index global_variable_index) const
{
  return get_relevant_feeval_vector<Rank>()[global_variable_index]
    .template get<type>()
    .get_hessian(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline DEAL_II_ALWAYS_INLINE
  typename FieldContainer<dim, degree, number>::template Hessian<Rank>
  FieldContainer<dim, degree, number>::get_hessian(Types::Index   global_variable_index,
                                                   DependencyType type) const
{
  return get_relevant_feeval_vector<Rank>()[global_variable_index].get(type).get_hessian(
    q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank, DependencyType type>
inline DEAL_II_ALWAYS_INLINE
  typename FieldContainer<dim, degree, number>::template Gradient<Rank>
  FieldContainer<dim, degree, number>::get_hessian_diagonal(
    Types::Index global_variable_index) const
{
  return get_relevant_feeval_vector<Rank>()[global_variable_index]
    .template get<type>()
    .get_hessian_diagonal(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline DEAL_II_ALWAYS_INLINE
  typename FieldContainer<dim, degree, number>::template Gradient<Rank>
  FieldContainer<dim, degree, number>::get_hessian_diagonal(
    Types::Index   global_variable_index,
    DependencyType type) const
{
  return get_relevant_feeval_vector<Rank>()[global_variable_index]
    .get(type)
    .get_hessian_diagonal(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank, DependencyType type>
inline DEAL_II_ALWAYS_INLINE
  typename FieldContainer<dim, degree, number>::template Value<Rank>
  FieldContainer<dim, degree, number>::get_laplacian(
    Types::Index global_variable_index) const
{
  return get_relevant_feeval_vector<Rank>()[global_variable_index]
    .template get<type>()
    .get_laplacian(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline DEAL_II_ALWAYS_INLINE
  typename FieldContainer<dim, degree, number>::template Value<Rank>
  FieldContainer<dim, degree, number>::get_laplacian(Types::Index   global_variable_index,
                                                     DependencyType type) const
{
  return get_relevant_feeval_vector<Rank>()[global_variable_index]
    .get(type)
    .get_laplacian(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank, DependencyType type>
inline DEAL_II_ALWAYS_INLINE typename FieldContainer<dim, degree, number>::ScalarValue
FieldContainer<dim, degree, number>::get_divergence(
  Types::Index global_variable_index) const
{
  static_assert(Rank == 1, "Divergences are only available for vector fields");
  return get_relevant_feeval_vector<Rank>()[global_variable_index]
    .template get<type>()
    .get_divergence(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline DEAL_II_ALWAYS_INLINE typename FieldContainer<dim, degree, number>::ScalarValue
FieldContainer<dim, degree, number>::get_divergence(Types::Index   global_variable_index,
                                                    DependencyType type) const
{
  static_assert(Rank == 1, "Divergences are only available for vector fields");
  return get_relevant_feeval_vector<Rank>()[global_variable_index]
    .get(type)
    .get_divergence(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank, DependencyType type>
inline DEAL_II_ALWAYS_INLINE dealii::
  SymmetricTensor<2, dim, typename FieldContainer<dim, degree, number>::ScalarValue>
  FieldContainer<dim, degree, number>::get_symmetric_gradient(
    Types::Index global_variable_index) const
{
  static_assert(Rank == 1, "Symmetric gradients are only available for vector fields");
  return get_relevant_feeval_vector<Rank>()[global_variable_index]
    .template get<type>()
    .get_symmetric_gradient(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline DEAL_II_ALWAYS_INLINE dealii::
  SymmetricTensor<2, dim, typename FieldContainer<dim, degree, number>::ScalarValue>
  FieldContainer<dim, degree, number>::get_symmetric_gradient(
    Types::Index   global_variable_index,
    DependencyType type) const
{
  static_assert(Rank == 1, "Symmetric gradients are only available for vector fields");
  return get_relevant_feeval_vector<Rank>()[global_variable_index]
    .get(type)
    .get_symmetric_gradient(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank, DependencyType type>
inline DEAL_II_ALWAYS_INLINE
  dealii::Tensor<1,
                 (dim == 2 ? 1 : dim),
                 typename FieldContainer<dim, degree, number>::ScalarValue>
  FieldContainer<dim, degree, number>::get_curl(Types::Index global_variable_index) const
{
  static_assert(Rank == 1, "Curl is only available for vector fields");
  Assert(dim > 1,
         dealii::ExcMessage(
           "Curl is only available for vector fields with dimension greater than 1."));
  return get_relevant_feeval_vector<Rank>()[global_variable_index]
    .template get<type>()
    .get_curl(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline DEAL_II_ALWAYS_INLINE
  dealii::Tensor<1,
                 (dim == 2 ? 1 : dim),
                 typename FieldContainer<dim, degree, number>::ScalarValue>
  FieldContainer<dim, degree, number>::get_curl(Types::Index   global_variable_index,
                                                DependencyType type) const
{
  static_assert(Rank == 1, "Curl is only available for vector fields");
  Assert(dim > 1,
         dealii::ExcMessage(
           "Curl is only available for vector fields with dimension greater than 1."));
  return get_relevant_feeval_vector<Rank>()[global_variable_index].get(type).get_curl(
    q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
inline DEAL_II_ALWAYS_INLINE
  dealii::Point<dim, typename FieldContainer<dim, degree, number>::ScalarValue>
  FieldContainer<dim, degree, number>::get_q_point_location() const
{
  return shared_feeval_scalar.quadrature_point(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
inline DEAL_II_ALWAYS_INLINE typename FieldContainer<dim, degree, number>::ScalarValue
FieldContainer<dim, degree, number>::get_element_volume() const
{
  return shared_feeval_scalar.JxW(q_point) /
         SystemWide<dim, degree>::quadrature.weight(q_point);
}

template <unsigned int dim, unsigned int degree, typename number>
inline DEAL_II_ALWAYS_INLINE unsigned int
FieldContainer<dim, degree, number>::get_n_q_points() const
{
  return shared_feeval_scalar.n_q_points;
}

template <unsigned int dim, unsigned int degree, typename number>
template <typename ValType>
inline DEAL_II_ALWAYS_INLINE void
FieldContainer<dim, degree, number>::set_value_term(Types::Index   global_variable_index,
                                                    const ValType &val)
{
  auto &relevant_feeval_vector =
    get_relevant_feeval_vector<RankFromVal<ValType>>()[global_variable_index];
  relevant_feeval_vector.template get<DependencyType::DST>().submit_value(val, q_point);
  relevant_feeval_vector.integration_flags |= dealii::EvaluationFlags::values;
}

template <unsigned int dim, unsigned int degree, typename number>
template <typename GradType>
inline DEAL_II_ALWAYS_INLINE void
FieldContainer<dim, degree, number>::set_gradient_term(Types::Index global_variable_index,
                                                       const GradType &val)
{
  auto &relevant_feeval_vector =
    get_relevant_feeval_vector<RankFromGrad<GradType>>()[global_variable_index];
  relevant_feeval_vector.template get<DependencyType::DST>().submit_gradient(val,
                                                                             q_point);
  relevant_feeval_vector.integration_flags |= dealii::EvaluationFlags::gradients;
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline DEAL_II_ALWAYS_INLINE std::vector<
  typename FieldContainer<dim, degree, number>::template FEEValuationDeps<Rank>> &
FieldContainer<dim, degree, number>::get_relevant_feeval_vector()
{
  if constexpr (Rank == TensorRank::Scalar)
    {
      return feeval_deps_scalar;
    }
  else if constexpr (Rank == TensorRank::Vector)
    {
      return feeval_deps_vector;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
template <TensorRank Rank>
inline DEAL_II_ALWAYS_INLINE const std::vector<
  typename FieldContainer<dim, degree, number>::template FEEValuationDeps<Rank>> &
FieldContainer<dim, degree, number>::get_relevant_feeval_vector() const
{
  if constexpr (Rank == TensorRank::Scalar)
    {
      return feeval_deps_scalar;
    }
  else if constexpr (Rank == TensorRank::Vector)
    {
      return feeval_deps_vector;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
inline DEAL_II_ALWAYS_INLINE void
FieldContainer<dim, degree, number>::feevaluation_size_valid(
  [[maybe_unused]] Types::Index field_index) const
{
  // TODO
}

template <unsigned int dim, unsigned int degree, typename number>
inline DEAL_II_ALWAYS_INLINE void
FieldContainer<dim, degree, number>::feevaluation_exists(
  [[maybe_unused]] Types::Index field_index,
  [[maybe_unused]] Types::Index dependency_index) const
{
  // TODO
}

template <unsigned int dim, unsigned int degree, typename number>
inline DEAL_II_ALWAYS_INLINE void
FieldContainer<dim, degree, number>::access_valid(
  [[maybe_unused]] Types::Index                             field_index,
  [[maybe_unused]] DependencyType                           dependency_type,
  [[maybe_unused]] dealii::EvaluationFlags::EvaluationFlags flag) const
{
  // TODO
}

template <unsigned int dim, unsigned int degree, typename number>
inline DEAL_II_ALWAYS_INLINE void
FieldContainer<dim, degree, number>::submission_valid(
  [[maybe_unused]] Types::Index   field_index,
  [[maybe_unused]] DependencyType dependency_type) const
{
  // TODO
}

PRISMS_PF_END_NAMESPACE
