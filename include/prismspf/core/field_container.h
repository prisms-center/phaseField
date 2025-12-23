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

#include <prismspf/core/solution_indexer.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include "prismspf/core/dependencies.h"
#include "prismspf/core/field_attributes.h"
#include "prismspf/core/solve_group.h"
#include "prismspf/core/variable_attributes.h"

#include <memory>
#include <type_traits>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class ElementVolume;

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
    static constexpr TensorRank rank_from_val  = TensorRank(Type::rank);
    static constexpr TensorRank rank_from_grad = TensorRank(Type::rank - 1);
  };

  template <>
  struct GetRankHelper<ScalarValue>
  {
    static constexpr TensorRank rank_from_val = TensorRank::Scalar;
    // TODO: make sure the enable_if is correct syntax
    template <typename std::enable_if<dim == 1>>
    static constexpr TensorRank rank_from_grad = TensorRank::Scalar;
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
  FieldContainer(
    const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
    const SolveGroup                                                       &_solve_group,
    const ElementVolume<dim, degree, number> &_element_volume,
    const std::vector<Types::Index>          &_global_to_local_solution,
    EquationType                              _equation_side);

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
  eval();

  /**
   * @brief Initialize based on cell, read solution vector, and evaluate based on
   * dependency flags for all dependencies.
   */
  void
  reinit_and_eval(unsigned int cell);

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
   * @brief Set the residual value of the specified scalar/vector field.
   */
  template <typename ValType, DependencyType type = DependencyType::Normal>
  void
  set_value_term(Types::Index global_variable_index, const ValType &val)
  {
    get_relevant_feeval_vector<GetRankFromVal<ValType>>()[global_variable_index]
      .template get<type>()
      .submit_value(val, q_point);
  }

  /**
   * @brief Set the residual gradient of the specified scalar/vector field.
   */
  template <typename GradType, DependencyType type = DependencyType::Normal>
  void
  set_gradient_term(Types::Index global_variable_index, const GradType &val)
  {
    get_relevant_feeval_vector<GetRankFromGrad<GradType>>()[global_variable_index]
      .template get<type>()
      .submit_gradient(val, q_point);
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
   * @brief Return the number of quadrature points.
   */
  [[nodiscard]] unsigned int
  get_n_q_points() const;

  /**
   * @brief Return the quadrature point location.
   */
  [[nodiscard]] dealii::Point<dim, ScalarValue>
  get_q_point_location() const;

  /**
   * @brief Integrate the residuals.
   */
  void
  integrate();

  /**
   * @brief Integrate the residuals and distribute from local to global.
   */
  void
  integrate_and_distribute();

  /**
   * @brief Struct to hold feevaluation relevant for this solve.
   */
  template <typename FEEvaluationType>
  struct FEEValuationDeps
  {
    using FEEDepPairPtr = std::unique_ptr<std::pair<FEEvaluationType, EvalFlags>>;
    FEEDepPairPtr                                            fe_eval;
    FEEDepPairPtr                                            fe_eval_change;
    std::array<FEEDepPairPtr, Numbers::max_saved_increments> fe_eval_old;
    /**
     * @brief The solution group and the block index
     * @note It would look nicer to just use the SolutionIndexer, but this way decreases
     * indexing
     */
    const SolutionLevel<dim, number> *solution_level = nullptr;
    unsigned int                      block_index    = -1;

    FEEValuationDeps(
      const Dependencies::Dependency                             &dependency,
      std::pair<const SolutionLevel<dim, number> *, unsigned int> mf_id_pair)
      : solution_level(mf_id_pair.first)
      , block_index(mf_id_pair.second)
    {
      if (dependency.flag)
        {
          fe_eval =
            std::make_unique<FEEvaluationType>(solution_level->matrix_free, block_index);
        }
      if (dependency.change_flag)
        {
          fe_eval_change =
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
    }

    template <DependencyType type>
    FEEvaluationType &
    get()
    {
      if constexpr (type == DependencyType::Normal)
        {
          return *fe_eval;
        }
      else if constexpr (type == DependencyType::Change)
        {
          return *fe_eval_change;
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
      if (fe_eval_change)
        {
          fe_eval_change->first.reinit(cell);
        }
      for (auto &old_fe_eval : fe_eval_old)
        {
          if (old_fe_eval)
            {
              old_fe_eval->first.reinit(cell);
            }
        }
    }

    void
    eval()
    {
      if (fe_eval)
        {
          fe_eval->first.read_dof_values_plain(
            solution_level->solutions.block(block_index));
          fe_eval->first.evaluate(fe_eval->second);
        }
      if (fe_eval_change)
        {
          fe_eval_change->first.read_dof_values_plain(
            solution_level->change_solutions.block(block_index));
          fe_eval_change->first.evaluate(fe_eval_change->second);
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
    }

    void
    reinit_and_eval(unsigned int cell)
    {
      if (fe_eval)
        {
          fe_eval->first.reinit(cell);
          fe_eval->first.read_dof_values_plain(
            solution_level->solutions.block(block_index));
          fe_eval->first.evaluate(fe_eval->second);
        }
      if (fe_eval_change)
        {
          fe_eval_change->first.reinit(cell);
          fe_eval_change->first.read_dof_values_plain(
            solution_level->change_solutions.block(block_index));
          fe_eval_change->first.evaluate(fe_eval_change->second);
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
    }
  };

  const std::vector<FieldAttributes>               *field_attributes_ptr;
  const SolveGroup                                 *solve_group;
  const SolutionIndexer<dim, number>               *solution_indexer;
  unsigned int                                      relative_level;
  std::vector<FEEValuationDeps<ScalarFEEvaluation>> feeval_deps_scalar;
  std::vector<FEEValuationDeps<VectorFEEvaluation>> feeval_deps_vector;
  std::vector<std::unique_ptr<ScalarFEEvaluation>>  dst_feeval_scalar;
  std::vector<std::unique_ptr<VectorFEEvaluation>>  dst_feeval_vector;

  /**
   * @brief The element volume container
   */
  const ElementVolume<dim, degree, number> *element_volume_handler;

  /**
   * @brief The quadrature point index.
   */
  unsigned int q_point = 0;

  /**
   * @brief Number of DoFs per cell.
   */
  unsigned int n_dofs_per_cell = 0;
};

PRISMS_PF_END_NAMESPACE
