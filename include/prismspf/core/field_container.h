// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

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
#include <prismspf/core/solve_group.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <memory>
#include <type_traits>
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
 */
template <unsigned int dim, unsigned int degree, typename number>
class FieldContainer
{
public:
  using TensorRank = FieldInfo::TensorRank;
  /**
   * @brief Typedef for the basic vector that apply our operations to.
   */
  using BlockVector    = SolutionIndexer<dim, number>::BlockVector;
  using SolutionVector = SolutionIndexer<dim, number>::SolutionVector;
  using MatrixFree     = dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>;

  /**
   * @brief Typedef for the basic value that the use manipulates.
   */
  using ScalarValue    = dealii::VectorizedArray<number>;
  using VectorValue    = dealii::Tensor<1, dim, ScalarValue>;
  using ScalarGradient = VectorValue;
  using VectorGradient = dealii::Tensor<2, dim, ScalarValue>;
  /**
   * @brief Typedef for scalar evaluation objects.
   */
  using ScalarFEEvaluation = dealii::
    FEEvaluation<dim, degree, degree + 1, 1, number, dealii::VectorizedArray<number>>;
  /**
   * @brief Typedef for vector evaluation objects.
   */
  using VectorFEEvaluation = dealii::
    FEEvaluation<dim, degree, degree + 1, dim, number, dealii::VectorizedArray<number>>;

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

  template <typename Type>
  struct GetRankHelper
  {
    static constexpr TensorRank rank_from_val = []() constexpr
      {
        if constexpr (std::is_same_v<Type, ScalarValue>)
          {
            return TensorRank::Scalar;
          }
        else
          {
            return TensorRank(Type::rank);
          }
      }();

    static constexpr TensorRank rank_from_grad = []() constexpr
      {
        if constexpr (std::is_same_v<Type, ScalarValue>) //&& dim == 1)
          {
            return TensorRank::Scalar;
          }
        else
          {
            return TensorRank(Type::rank - 1);
          }
      }();
  };

  template <typename ValType>
  using GetRankFromVal = GetRankHelper<ValType>::rank_from_val;
  template <typename GradType>
  using GetRankFromGrad = GetRankHelper<GradType>::rank_from_grad;

  template <TensorRank Rank>
  std::vector<FEEval<Rank>>
  get_relevant_feeval_vector()
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

  /**
   * @brief Typedef for scalar diagonal matrix objects.
   */
  using ScalarDiagonal = dealii::AlignedVector<ScalarValue>;

  /**
   * @brief Typedef for vector diagonal matrix objects.
   */
  using VectorDiagonal = dealii::AlignedVector<VectorValue>;

  /**
   * @brief Constructor.
   */
  FieldContainer(const std::vector<FieldAttributes> &_field_attributes,
                 const SolutionIndexer<dim, number> &_solution_indexer,
                 unsigned int                        _relative_level,
                 const DependencySet                &dependency_map,
                 const SolveGroup                   &_solve_group,
                 const MatrixFree                   &matrix_free);

  /**
   * @brief Initialize based on cell for all dependencies.
   */
  void
  reinit(unsigned int cell);

  /**
   * @brief Read solution vector, and evaluate based on
   * dependency flags for all dependencies.
   */
  void
  eval(const BlockVector *src_solutions);

  /**
   * @brief Initialize based on cell, read solution vector, and evaluate based on
   * dependency flags for all dependencies.
   */
  void
  reinit_and_eval(unsigned int cell, const BlockVector *src_solutions);

  /**
   * @brief Integrate the residuals.
   */
  void
  integrate();

  /**
   * @brief Integrate the residuals and distribute from local to global.
   */
  void
  integrate_and_distribute(BlockVector *dst_solutions);

  /**
   * @brief Set the current quadtrature point
   */
  void
  set_q_point(unsigned int quad)
  {
    q_point = quad;
  }

  /**
   * @brief Return the value of the specified field.
   *
   * @tparam T the return type. Must be either a `ScalarValue` or `dealii::Tensor<1, dim,
   * ScalarValue>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <TensorRank Rank, DependencyType type = DependencyType::Normal>
  [[nodiscard]] Value<Rank>
  get_value(Types::Index global_variable_index) const
  {
    return get_relevant_feeval_vector<Rank>()[global_variable_index]
      .template get<type>()
      .get_value(q_point);
  }

  /**
   * @brief Return the gradient of the specified field.
   *
   * @tparam T the return type. Must be either a `dealii::Tensor<1, dim,
   * ScalarValue>` or `dealii::Tensor<2, dim, ScalarValue>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <TensorRank Rank, DependencyType type = DependencyType::Normal>
  [[nodiscard]] Gradient<Rank>
  get_gradient(Types::Index global_variable_index) const
  {
    return get_relevant_feeval_vector<Rank>()[global_variable_index]
      .template get<type>()
      .get_gradient(q_point);
  }

  /**
   * @brief Return the hessian of the specified field.
   *
   * @tparam T the return type. Must be either a `dealii::Tensor<2, dim,
   * ScalarValue>` or `dealii::Tensor<3, dim, ScalarValue>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <TensorRank Rank, DependencyType type = DependencyType::Normal>
  [[nodiscard]] Hessian<Rank>
  get_hessian(Types::Index global_variable_index) const
  {
    return get_relevant_feeval_vector<Rank>()[global_variable_index]
      .template get<type>()
      .get_hessian(q_point);
  }

  /**
   * @brief Return the diagonal of the hessian of the specified field.
   *
   * @tparam T the return type. Must be either a `dealii::Tensor<1, dim,
   * ScalarValue>` or `dealii::Tensor<2, dim, ScalarValue>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <TensorRank Rank, DependencyType type = DependencyType::Normal>
  [[nodiscard]] Gradient<Rank>
  get_hessian_diagonal(Types::Index global_variable_index) const
  {
    return get_relevant_feeval_vector<Rank>()[global_variable_index]
      .template get<type>()
      .get_hessian_diagonal(q_point);
  }

  /**
   * @brief Return the laplacian of the specified field.
   *
   * @tparam T the return type. Must be either a `ScalarValue` or `dealii::Tensor<1, dim,
   * ScalarValue>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <TensorRank Rank, DependencyType type = DependencyType::Normal>
  [[nodiscard]] Value<Rank>
  get_laplacian(Types::Index global_variable_index) const
  {
    return get_relevant_feeval_vector<Rank>()[global_variable_index]
      .template get<type>()
      .get_laplacian(q_point);
  }

  /**
   * @brief Return the divergence of the specified field.
   *
   * @tparam T the return type. Must be `ScalarValue`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  // TODO: FIgure out the assertion here. Dealii probly will have one, but we should too.
  template <TensorRank Rank, DependencyType type = DependencyType::Normal>
  [[nodiscard]] auto
  get_divergence(Types::Index global_variable_index) const
    -> ScalarValue /* Value<Rank>::value_type */
  {
    return get_relevant_feeval_vector<Rank>()[global_variable_index]
      .template get<type>()
      .get_divergence(q_point);
  }

  /**
   * @brief Return the symmetric gradient of the specified field.
   *
   * @tparam T the return type. Must be either a `dealii::Tensor<2, dim,
   * ScalarValue>` or `dealii::SymmetricTensor<2, dim, ScalarValue>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  // TODO: FIgure out the assertion here. Dealii probly will have one, but we should too.
  template <TensorRank Rank, DependencyType type = DependencyType::Normal>
  [[nodiscard]] auto
  get_symmetric_gradient(Types::Index global_variable_index) const
    -> dealii::SymmetricTensor<2, dim, ScalarValue>
  {
    return get_relevant_feeval_vector<Rank>()[global_variable_index]
      .template get<type>()
      .get_symmetric_gradient(q_point);
  }

  /**
   * @brief Return the curl of the specified field.
   *
   * Note that this is dealii::VectorizedArray<number> type for 2D and dealii::Tensor<1,
   * dim, dealii::VectorizedArray<number>> type for 3D.
   *
   * @tparam T the return type. Must be either a `dealii::Tensor<1, 1,
   * ScalarValue>` or `dealii::Tensor<1, dim, ScalarValue>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  // TODO: FIgure out the assertion here. Dealii probly will have one, but we should too.
  template <TensorRank Rank, DependencyType type = DependencyType::Normal>
  [[nodiscard]] auto
  get_curl(Types::Index global_variable_index) const
    -> dealii::Tensor<1, (dim == 2 ? 1 : dim), ScalarValue>
  {
    return get_relevant_feeval_vector<Rank>()[global_variable_index]
      .template get<type>()
      .get_curl(q_point);
  }

  /**
   * @brief Return the quadrature point location.
   */
  [[nodiscard]] dealii::Point<dim, ScalarValue>
  get_q_point_location() const
  {
    return shared_feeval_scalar.quadrature_point(q_point);
  }

  /**
   * @brief Return the number of quadrature points.
   */
  [[nodiscard]] unsigned int
  get_n_q_points() const
  {
    return shared_feeval_scalar.n_q_points;
  }

  /**
   * @brief Set the residual value of the specified scalar/vector field.
   */
  template <typename ValType>
  void
  set_value_term(Types::Index global_variable_index, const ValType &val)
  {
    auto &relevant_feeval_vector =
      get_relevant_feeval_vector<GetRankFromVal<ValType>>()[global_variable_index];
    relevant_feeval_vector.template get<DependencyType::Change>().submit_value(val,
                                                                               q_point);
    relevant_feeval_vector.integration_flags |= dealii::EvaluationFlags::values;
  }

  /**
   * @brief Set the residual gradient of the specified scalar/vector field.
   */
  template <typename GradType>
  void
  set_gradient_term(Types::Index global_variable_index, const GradType &val)
  {
    auto &relevant_feeval_vector =
      get_relevant_feeval_vector<GetRankFromGrad<GradType>>()[global_variable_index];
    relevant_feeval_vector.template get<DependencyType::Change>()
      .submit_gradient(val, q_point);
    relevant_feeval_vector.integration_flags |= dealii::EvaluationFlags::gradients;
  }

private:
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
   * @brief Struct to hold feevaluation relevant for this solve.
   */
  template <typename FEEvaluationType>
  struct FEEValuationDeps
  {
    using FEEDepPairPtr = std::unique_ptr<std::pair<FEEvaluationType, EvalFlags>>;
    FEEDepPairPtr                                            fe_eval;
    FEEDepPairPtr                                            fe_eval_src_dst;
    std::array<FEEDepPairPtr, Numbers::max_saved_increments> fe_eval_old;
    EvalFlags integration_flags = EvalFlags::nothing;
    /**
     * @brief The solution group and the block index
     * @note It would look nicer to just use the SolutionIndexer, but this way decreases
     * indexing
     */
    const SolutionLevel<dim, number> *solution_level = nullptr;
    unsigned int                      block_index    = -1;

    FEEValuationDeps(
      const Dependencies::Dependency                             &dependency,
      std::pair<const SolutionLevel<dim, number> *, unsigned int> mf_id_pair,
      bool                                                        is_dst)
      : solution_level(mf_id_pair.first)
      , block_index(mf_id_pair.second)
    {
      if (dependency.flag)
        {
          fe_eval =
            std::make_unique<FEEvaluationType>(solution_level->matrix_free, block_index);
        }
      for (unsigned int age = 0; age < Numbers::max_saved_increments; ++age)
        {
          if (dependency.old_flags.at(age))
            {
              fe_eval_old[age] =
                std::make_unique<FEEvaluationType>(solution_level->matrix_free,
                                                   block_index);
            }
        }
      if (dependency.change_flag || is_dst)
        {
          fe_eval_src_dst =
            std::make_unique<FEEvaluationType>(solution_level->matrix_free, block_index);
        }
    }

    template <DependencyType type>
    FEEvaluationType &
    get()
    {
      // TODO: Assertions
      if constexpr (type == DependencyType::Change)
        {
          return *fe_eval_src_dst;
        }
      else if constexpr (type == DependencyType::Normal)
        {
          return *fe_eval;
        }
      else
        {
          return *fe_eval_old[int(type)];
        }
    }

    void
    reinit(unsigned int cell)
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

    void
    eval(const BlockVector *_src_solutions)
    {
      if (fe_eval)
        {
          fe_eval->first.read_dof_values_plain(
            solution_level->solutions.block(block_index));
          fe_eval->first.evaluate(fe_eval->second);
        }
      for (unsigned int age = 0; age < Numbers::max_saved_increments; ++age)
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
          fe_eval_src_dst->first.read_dof_values_plain(
            _src_solutions->block(block_index));
          fe_eval_src_dst->first.evaluate(fe_eval_src_dst->second);
        }
    }

    void
    reinit_and_eval(unsigned int cell, const BlockVector *_src_solutions)
    {
      if (fe_eval)
        {
          fe_eval->first.reinit(cell);
          fe_eval->first.read_dof_values_plain(
            solution_level->solutions.block(block_index));
          fe_eval->first.evaluate(fe_eval->second);
        }
      for (unsigned int age = 0; age < Numbers::max_saved_increments; ++age)
        {
          if (FEEDepPairPtr &old_fe_eval = fe_eval_old[age])
            {
              old_fe_eval->first.reinit(cell);
              old_fe_eval->first.read_dof_values_plain(
                solution_level->old_solutions[age].block(block_index));
              old_fe_eval->first.evaluate(old_fe_eval->second);
            }
        }
      if (fe_eval_src_dst && fe_eval_src_dst->second != EvalFlags::nothing)
        {
          fe_eval_src_dst->first.reinit(cell);
          fe_eval_src_dst->first.read_dof_values_plain(
            _src_solutions->block(block_index));
          fe_eval_src_dst->first.evaluate(fe_eval_src_dst->second);
        }
    }

    void
    integrate()
    {
      if (fe_eval_src_dst)
        {
          fe_eval_src_dst->first.integrate(integration_flags);
        }
    }

    void
    integrate_and_distribute(BlockVector &dst_solutions)
    {
      if (fe_eval_src_dst)
        {
          fe_eval_src_dst->first.integrate_scatter(integration_flags,
                                                   dst_solutions.block(block_index));
        }
    }
  };

  //================================================================================
  // Members relevant for accessing the all fields and dependencies
  //================================================================================
  const std::vector<FieldAttributes>               *field_attributes_ptr;
  const SolutionIndexer<dim, number>               *solution_indexer;
  std::vector<FEEValuationDeps<ScalarFEEvaluation>> feeval_deps_scalar;
  std::vector<FEEValuationDeps<VectorFEEvaluation>> feeval_deps_vector;

  //================================================================================
  // Members relevant for submitting values and integrating and accessing lhs fields
  //================================================================================
  const SolveGroup *solve_group;

  //================================================================================
  // Shared members for generic access
  //================================================================================
  ScalarFEEvaluation shared_feeval_scalar;
  unsigned int       relative_level;

  /**
   * @brief The quadrature point index.
   */
  unsigned int q_point = 0;
};

PRISMS_PF_END_NAMESPACE
