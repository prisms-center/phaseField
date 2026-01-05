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
class DSTContainer
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
        return feeval_scalar;
      }
    else if constexpr (Rank == TensorRank::Vector)
      {
        return feeval_vector;
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
  DSTContainer(const dealii::MatrixFree<dim, number, ScalarValue> &data,
               const SolveGroup                                   &_solve_group,
               const std::vector<FieldAttributes>                 &_field_attributes,
               SolutionIndexer<dim, number>                       &_solution_indexer);

  /**
   * @brief Initialize based on cell for all dependencies.
   */
  void
  reinit(unsigned int cell);

  /**
   * @brief Set the residual value of the specified scalar/vector field.
   */
  template <typename ValType, DependencyType type = DependencyType::Normal>
  void
  set_value_term(Types::Index global_variable_index, const ValType &val)
  {
    get_relevant_feeval_vector<GetRankFromVal<ValType>>()[global_variable_index]
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
      .submit_gradient(val, q_point);
  }

  /**
   * @brief Check that a value is valid for submission.
   */
  void
  submission_valid(Types::Index field_index, DependencyType dependency_type) const;

private:
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
    FEEDepPairPtr fe_eval;
    /**
     * @brief The solution group and the block index
     * @note It would look nicer to just use the SolutionIndexer, but this way decreases
     * indexing
     */
    const SolutionLevel<dim, number> *solution_level = nullptr;
    unsigned int                      block_index    = -1;
  };

  const std::vector<FieldAttributes>              *field_attributes_ptr;
  const SolveGroup                                *solve_group;
  const SolutionIndexer<dim, number>              *solution_indexer;
  unsigned int                                     relative_level;
  std::vector<std::unique_ptr<ScalarFEEvaluation>> feeval_scalar;
  std::vector<std::unique_ptr<VectorFEEvaluation>> feeval_vector;

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
