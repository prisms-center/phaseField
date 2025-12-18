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
#include "prismspf/core/solve_group.h"
#include "prismspf/core/variable_attributes.h"

#include <memory>
#include <type_traits>
#include <variant>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

struct VariableAttributes;

template <unsigned int dim, unsigned int degree, typename number>
class ElementVolume;

// clang-format off

/**
 * @brief Overload pattern for lambdas.
 */
template <class... Ts>
struct Overload : Ts...{ using Ts::operator()...; };
template <class... Ts>
Overload(Ts...) -> Overload<Ts...>;

// clang-format on

enum class EquationType : int
{
  RHS,
  LHS,
  Residual = RHS
};

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
    template <typename = std::enable_if<dim == 1>>
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
   * @brief Typedef for the variant evaluation objects. Note that the states become
   * degenerate at dim = 1, hence the use of std::conditional_t.
   */
  using VariantFEEvaluation =
    std::conditional_t<std::is_same_v<ScalarFEEvaluation, VectorFEEvaluation>,
                       std::unique_ptr<ScalarFEEvaluation>,
                       std::variant<std::unique_ptr<ScalarFEEvaluation>,
                                    std::unique_ptr<VectorFEEvaluation>>>;

  /**
   * @brief Typedef for scalar diagonal matrix objects.
   */
  using ScalarDiagonal = dealii::AlignedVector<ScalarValue>;

  /**
   * @brief Typedef for vector diagonal matrix objects.
   */
  using VectorDiagonal = dealii::AlignedVector<VectorValue>;

  /**
   * @brief Typedef for the variant diagonal matrix objects.
   */
  using VariantDiagonal =
    std::variant<std::unique_ptr<ScalarDiagonal>, std::unique_ptr<VectorDiagonal>>;
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
  feevaluation_size_valid(Types::Index field_index, Types::Index dependency_index) const;

  /**
   * @brief Check whether the entry for the FEEvaluation is within the bounds of the
   * vector and not a nullptr.
   */
  void
  feevaluation_exists(Types::Index field_index, Types::Index dependency_index) const;

  /**
   * @brief Check whether the entry for the global solution vector to local one is within
   * the bounds of the vector and contains a valid entry.
   */
  void
  global_to_local_solution_exists(Types::Index field_index,
                                  Types::Index dependency_index) const;

  /**
   * @brief Get the local solution index for the given field index and dependency type.
   */
  [[nodiscard]] Types::Index
  get_local_solution_index(Types::Index field_index, Types::Index dependency_index) const;

  /**
   * @brief Get the local solution index for the given field index and dependency type.
   */
  [[nodiscard]] Types::Index
  get_local_solution_index(Types::Index   field_index,
                           DependencyType dependency_type) const;

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
  submission_valid(DependencyType dependency_type) const;

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
   * @brief Initialize, read DOFs, and set evaulation flags for each variable.
   */
  void
  reinit_and_eval(const std::vector<SolutionVector *> &src, unsigned int cell);

  /**
   * @brief Initialize, read DOFs, and set evaulation flags for each variable.
   */
  void
  reinit_and_eval(const SolutionVector &src, unsigned int cell);

  /**
   * @brief Initialize the cell for all dependencies of a certain variable index.
   */
  void
  reinit(unsigned int cell, Types::Index global_variable_index);

  /**
   * @brief Read dofs values on the cell for all dependencies of a certain variable index.
   */
  void
  read_dof_values(const std::vector<SolutionVector *> &src);

  /**
   * @brief Evaluate the flags on the cell for all dependencies of a certain variable
   * index.
   */
  void
  eval(Types::Index global_variable_index);

  /**
   * @brief Integrate the residuals and distribute from local to global.
   */
  void
  integrate_and_distribute(std::vector<SolutionVector *> &dst);

  /**
   * @brief Integrate the residuals and distribute from local to global.
   */
  void
  integrate_and_distribute(SolutionVector &dst);

  /**
   * @brief Integrate the residuals for a certain variable index.
   */
  void
  integrate(Types::Index global_variable_index);

  /**
   * @brief Get the FEEvaluation pointer from the variant.
   */
  template <typename FEEvaluationType>
  FEEvaluationType *
  extract_feeval_ptr(VariantFEEvaluation &variant) const;

  /**
   * @brief Get the diagonal pointer from the variant.
   */
  template <typename DiagonalType>
  DiagonalType *
  extract_diagonal_ptr(VariantDiagonal &variant) const;

  /**
   * @brief Evaluate the diagonal entry for a given cell.
   */
  template <typename FEEvaluationType, typename DiagonalType>
  void
  eval_cell_diagonal(FEEvaluationType *feeval_ptr,
                     DiagonalType     *diagonal_ptr,
                     unsigned int      cell,
                     Types::Index      global_variable_index,
                     const std::function<void(FieldContainer &,
                                              const dealii::Point<dim, ScalarValue> &,
                                              const ScalarValue &)> &func,
                     SolutionVector                                 &dst,
                     const std::vector<SolutionVector *>            &src_subset);

  template <typename FEEvaluationType>
  struct FEEValuationDeps
  {
    // Apparently FEEvaluation objects are cheap, so maybe it's overkill to use dynamic
    // mem here.
    std::unique_ptr<FEEvaluationType>                           fe_eval;
    std::unique_ptr<FEEvaluationType>                           fe_eval_change;
    std::array<FEEvaluationType, Numbers::max_saved_increments> fe_eval_old;

    FEEValuationDeps(const Dependencies::Dependency     &dependency,
                     std::pair<MatrixFree, unsigned int> mf_id_pair)
    {
      if (dependency.flag)
        {
          fe_eval =
            std::make_unique<FEEvaluationType>(mf_id_pair.first, mf_id_pair.second);
        }
      if (dependency.change_flag)
        {
          fe_eval_change =
            std::make_unique<FEEvaluationType>(mf_id_pair.first, mf_id_pair.second);
        }
      for (unsigned int age = 0; age < Numbers::max_saved_increments; ++age)
        {
          if (dependency.old_flags.at(age))
            {
              fe_eval_old[age] =
                std::make_unique<FEEvaluationType>(mf_id_pair.first, mf_id_pair.second);
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
  };

  const SolveGroup                                 *solve_group;
  EquationType                                      equation_type;
  const SolutionIndexer<dim, number>               *solution_indexer;
  std::vector<FEEValuationDeps<ScalarFEEvaluation>> feeval_deps_scalar;
  std::vector<FEEValuationDeps<VectorFEEvaluation>> feeval_deps_vector;
  unsigned int                                      relative_level;

  /**
   * @brief Max number of fields.
   */
  Types::Index max_fields = Numbers::invalid_index;

  /**
   * @brief Nax number of dependency types.
   */
  Types::Index max_dependency_types = Numbers::invalid_index;

  /**
   * @brief Vector of FEEvaluation objects for each active variable.
   *
   * The value is a variant that can hold either a ptr to a scalar or vector FEEvaluation.
   * For performance reasons, we have a vector with length of max_fields *
   * max_dependency_types. Consequently, most of the vector is filled with nullptr's.
   */
  std::vector<VariantFEEvaluation> feeval_vector;

  /**
   * @brief The attribute list of the relevant subset of variables.
   */
  const std::map<Types::Index, VariableAttributes> *subset_attributes;

  /**
   * @brief The element volume container
   */
  const ElementVolume<dim, degree, number> *element_volume_handler;

  /**
   * @brief Mapping from global solution vectors to the local ones
   */
  const std::vector<Types::Index> *global_to_local_solution;

  /**
   * @brief The solve type
   */
  SolveType solve_type;

  /**
   * @brief The quadrature point index.
   */
  unsigned int q_point = 0;

  /**
   * @brief Number of DoFs per cell.
   */
  unsigned int n_dofs_per_cell = 0;

  /**
   * @brief Diagonal matrix that is used for preconditioning of fields.
   */
  VariantDiagonal diagonal;
};

PRISMS_PF_END_NAMESPACE
