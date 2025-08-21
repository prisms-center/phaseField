// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <type_traits>
#include <variant>

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
class VariableContainer
{
public:
  /**
   * @brief Typedef for the basic vector that apply our operations to.
   */
  using VectorType = dealii::LinearAlgebra::distributed::Vector<number>;

  /**
   * @brief Typedef for the basic value that the use manipulates.
   */
  using SizeType = dealii::VectorizedArray<number>;

  /**
   * @brief Constructor.
   */
  VariableContainer(
    const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
    const std::map<Types::Index, VariableAttributes> &_subset_attributes,
    const ElementVolume<dim, degree, number>         &_element_volume,
    const std::vector<Types::Index>                  &_global_to_local_solution,
    const SolveType                                  &_solve_type,
    bool                                              use_local_mapping = false);

  /**
   * @brief Return the value of the specified field.
   *
   * @tparam T the return type. Must be either a `SizeType` or `dealii::Tensor<1, dim,
   * SizeType>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <typename T>
  [[nodiscard]] T
  get_value(Types::Index   global_variable_index,
            DependencyType dependency_type = DependencyType::Normal) const
  requires(std::is_same_v<T, SizeType> ||
           std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::values);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    // For the 1D case the FEEvaluation objects for scalar/vector fields are degenerate.
    // Otherwise, use std::visit
    if constexpr (dim == 1)
      {
        // Wrap the vector value here
        if constexpr (std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
          {
            dealii::Tensor<1, dim, SizeType> wrapper;
            wrapper[0] = feeval_vector[(global_variable_index * max_dependency_types) +
                                       static_cast<Types::Index>(dependency_type)]
                           ->get_value(q_point);
            return wrapper;
          }
        else
          {
            return feeval_vector[(global_variable_index * max_dependency_types) +
                                 static_cast<Types::Index>(dependency_type)]
              ->get_value(q_point);
          }
      }
    else
      {
        return std::visit<T>(
          [&](auto &feeval_ptr) -> T
          {
            using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;

            if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
              {
                if constexpr (std::is_same_v<T, SizeType>)
                  {
                    return feeval_ptr->get_value(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match scalar FEEvaluation"));
                    return T {};
                  }
              }
            else if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
              {
                if constexpr (std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
                  {
                    return feeval_ptr->get_value(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match vector FEEvaluation"));
                    return T {};
                  }
              }
          },
          feeval_vector[(global_variable_index * max_dependency_types) +
                        static_cast<Types::Index>(dependency_type)]);
      }
  }

  /**
   * @brief Return the gradient of the specified field.
   *
   * @tparam T the return type. Must be either a `dealii::Tensor<1, dim,
   * SizeType>` or `dealii::Tensor<2, dim, SizeType>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <typename T>
  [[nodiscard]] T
  get_gradient(Types::Index   global_variable_index,
               DependencyType dependency_type = DependencyType::Normal) const
  requires(std::is_same_v<T, dealii::Tensor<1, dim, SizeType>> ||
           std::is_same_v<T, dealii::Tensor<2, dim, SizeType>>)
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::gradients);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    // For the 1D case the FEEvaluation objects for scalar/vector fields are degenerate.
    // Otherwise, use std::visit
    if constexpr (dim == 1)
      {
        // Wrap the vector value here
        if constexpr (std::is_same_v<T, dealii::Tensor<2, dim, SizeType>>)
          {
            dealii::Tensor<2, dim, SizeType> wrapper;
            wrapper[0] = feeval_vector[(global_variable_index * max_dependency_types) +
                                       static_cast<Types::Index>(dependency_type)]
                           ->get_gradient(q_point);
            return wrapper;
          }
        else
          {
            return feeval_vector[(global_variable_index * max_dependency_types) +
                                 static_cast<Types::Index>(dependency_type)]
              ->get_gradient(q_point);
          }
      }
    else
      {
        return std::visit<T>(
          [&](auto &feeval_ptr) -> T
          {
            using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;

            if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
              {
                if constexpr (std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
                  {
                    return feeval_ptr->get_gradient(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match scalar FEEvaluation"));
                    return T {};
                  }
              }
            else if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
              {
                if constexpr (std::is_same_v<T, dealii::Tensor<2, dim, SizeType>>)
                  {
                    return feeval_ptr->get_gradient(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match vector FEEvaluation"));
                    return T {};
                  }
              }
          },
          feeval_vector[(global_variable_index * max_dependency_types) +
                        static_cast<Types::Index>(dependency_type)]);
      }
  }

  /**
   * @brief Return the hessian of the specified field.
   *
   * @tparam T the return type. Must be either a `dealii::Tensor<2, dim,
   * SizeType>` or `dealii::Tensor<3, dim, SizeType>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <typename T>
  [[nodiscard]] T
  get_hessian(Types::Index   global_variable_index,
              DependencyType dependency_type = DependencyType::Normal) const
  requires(std::is_same_v<T, dealii::Tensor<2, dim, SizeType>> ||
           std::is_same_v<T, dealii::Tensor<3, dim, SizeType>>)
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::hessians);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    // For the 1D case the FEEvaluation objects for scalar/vector fields are degenerate.
    // Otherwise, use std::visit
    if constexpr (dim == 1)
      {
        // Wrap the vector value here
        if constexpr (std::is_same_v<T, dealii::Tensor<3, dim, SizeType>>)
          {
            dealii::Tensor<3, dim, SizeType> wrapper;
            wrapper[0] = feeval_vector[(global_variable_index * max_dependency_types) +
                                       static_cast<Types::Index>(dependency_type)]
                           ->get_hessian(q_point);
            return wrapper;
          }
        else
          {
            return feeval_vector[(global_variable_index * max_dependency_types) +
                                 static_cast<Types::Index>(dependency_type)]
              ->get_hessian(q_point);
          }
      }
    else
      {
        return std::visit<T>(
          [&](auto &feeval_ptr) -> T
          {
            using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;

            if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
              {
                if constexpr (std::is_same_v<T, dealii::Tensor<2, dim, SizeType>>)
                  {
                    return feeval_ptr->get_hessian(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match scalar FEEvaluation"));
                    return T {};
                  }
              }
            else if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
              {
                if constexpr (std::is_same_v<T, dealii::Tensor<3, dim, SizeType>>)
                  {
                    return feeval_ptr->get_hessian(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match vector FEEvaluation"));
                    return T {};
                  }
              }
          },
          feeval_vector[(global_variable_index * max_dependency_types) +
                        static_cast<Types::Index>(dependency_type)]);
      }
  }

  /**
   * @brief Return the diagonal of the hessian of the specified field.
   *
   * @tparam T the return type. Must be either a `dealii::Tensor<1, dim,
   * SizeType>` or `dealii::Tensor<2, dim, SizeType>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <typename T>
  [[nodiscard]] T
  get_hessian_diagonal(Types::Index   global_variable_index,
                       DependencyType dependency_type = DependencyType::Normal) const
  requires(std::is_same_v<T, dealii::Tensor<1, dim, SizeType>> ||
           std::is_same_v<T, dealii::Tensor<2, dim, SizeType>>)
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::hessians);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    // For the 1D case the FEEvaluation objects for scalar/vector fields are degenerate.
    // Otherwise, use std::visit
    if constexpr (dim == 1)
      {
        // Wrap the vector value here
        if constexpr (std::is_same_v<T, dealii::Tensor<2, dim, SizeType>>)
          {
            dealii::Tensor<2, dim, SizeType> wrapper;
            wrapper[0] = feeval_vector[(global_variable_index * max_dependency_types) +
                                       static_cast<Types::Index>(dependency_type)]
                           ->get_hessian_diagonal(q_point);
            return wrapper;
          }
        else
          {
            return feeval_vector[(global_variable_index * max_dependency_types) +
                                 static_cast<Types::Index>(dependency_type)]
              ->get_hessian_diagonal(q_point);
          }
      }
    else
      {
        return std::visit<T>(
          [&](auto &feeval_ptr) -> T
          {
            using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;

            if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
              {
                if constexpr (std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
                  {
                    return feeval_ptr->get_hessian_diagonal(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match scalar FEEvaluation"));
                    return T {};
                  }
              }
            else if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
              {
                if constexpr (std::is_same_v<T, dealii::Tensor<2, dim, SizeType>>)
                  {
                    return feeval_ptr->get_hessian_diagonal(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match vector FEEvaluation"));
                    return T {};
                  }
              }
          },
          feeval_vector[(global_variable_index * max_dependency_types) +
                        static_cast<Types::Index>(dependency_type)]);
      }
  }

  /**
   * @brief Return the laplacian of the specified field.
   *
   * @tparam T the return type. Must be either a `SizeType` or `dealii::Tensor<1, dim,
   * SizeType>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <typename T>
  [[nodiscard]] T
  get_laplacian(Types::Index   global_variable_index,
                DependencyType dependency_type = DependencyType::Normal) const
  requires(std::is_same_v<T, SizeType> ||
           std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::hessians);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    // For the 1D case the FEEvaluation objects for scalar/vector fields are degenerate.
    // Otherwise, use std::visit
    if constexpr (dim == 1)
      {
        // Wrap the vector value here
        if constexpr (std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
          {
            dealii::Tensor<1, dim, SizeType> wrapper;
            wrapper[0] = feeval_vector[(global_variable_index * max_dependency_types) +
                                       static_cast<Types::Index>(dependency_type)]
                           ->get_laplacian(q_point);
            return wrapper;
          }
        else
          {
            return feeval_vector[(global_variable_index * max_dependency_types) +
                                 static_cast<Types::Index>(dependency_type)]
              ->get_laplacian(q_point);
          }
      }
    else
      {
        return std::visit<T>(
          [&](auto &feeval_ptr) -> T
          {
            using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;

            if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
              {
                if constexpr (std::is_same_v<T, SizeType>)
                  {
                    return feeval_ptr->get_laplacian(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match scalar FEEvaluation"));
                    return T {};
                  }
              }
            else if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
              {
                if constexpr (std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
                  {
                    return feeval_ptr->get_laplacian(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match vector FEEvaluation"));
                    return T {};
                  }
              }
          },
          feeval_vector[(global_variable_index * max_dependency_types) +
                        static_cast<Types::Index>(dependency_type)]);
      }
  }

  /**
   * @brief Return the divergence of the specified field.
   *
   * @tparam T the return type. Must be `SizeType`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <typename T>
  [[nodiscard]] T
  get_divergence(Types::Index   global_variable_index,
                 DependencyType dependency_type = DependencyType::Normal) const
  requires(std::is_same_v<T, SizeType>)
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::gradients);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    // For the 1D case the FEEvaluation objects for scalar/vector fields are degenerate.
    // Otherwise, use std::visit
    if constexpr (dim == 1)
      {
        // Typically, we wrap the vector value for consistency with higher space
        // dimensions; however, the divergence of a vector field (n_components = dim),
        // always returns a SizeType.
        return feeval_vector[(global_variable_index * max_dependency_types) +
                             static_cast<Types::Index>(dependency_type)]
          ->get_divergence(q_point);
      }
    else
      {
        return std::visit<T>(
          [&](auto &feeval_ptr) -> T
          {
            using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;

            if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
              {
                Assert(false,
                       dealii::ExcMessage(
                         "You cannot take the divergence of scalar fields."));
                return T {};
              }
            else if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
              {
                if constexpr (std::is_same_v<T, SizeType>)
                  {
                    return feeval_ptr->get_divergence(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match vector FEEvaluation"));
                    return T {};
                  }
              }
          },
          feeval_vector[(global_variable_index * max_dependency_types) +
                        static_cast<Types::Index>(dependency_type)]);
      }
  }

  /**
   * @brief Return the symmetric gradient of the specified field.
   *
   * @tparam T the return type. Must be either a `dealii::Tensor<2, dim,
   * SizeType>` or `dealii::SymmetricTensor<2, dim, SizeType>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <typename T>
  [[nodiscard]] T
  get_symmetric_gradient(Types::Index   global_variable_index,
                         DependencyType dependency_type = DependencyType::Normal) const
  requires(std::is_same_v<T, dealii::Tensor<2, dim, SizeType>> ||
           std::is_same_v<T, dealii::SymmetricTensor<2, dim, SizeType>>)
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::gradients);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    // For the 1D case the FEEvaluation objects for scalar/vector fields are degenerate.
    // Otherwise, use std::visit
    if constexpr (dim == 1)
      {
        // Typically, we wrap the vector value for consistency with higher space
        // dimensions; however, the symmetric gradient of a vector field (n_components =
        // dim), always returns a dealii::SymmetricTensor<2, dim, SizeType>.
        return feeval_vector[(global_variable_index * max_dependency_types) +
                             static_cast<Types::Index>(dependency_type)]
          ->get_symmetric_gradient(q_point);
      }
    else
      {
        return std::visit<T>(
          [&](auto &feeval_ptr) -> T
          {
            using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;

            if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
              {
                Assert(false,
                       dealii::ExcMessage(
                         "You cannot take the symmetric gradient of scalar fields."));
                return T {};
              }
            else if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
              {
                if constexpr (std::is_same_v<T, dealii::Tensor<2, dim, SizeType>> ||
                              std::is_same_v<T,
                                             dealii::SymmetricTensor<2, dim, SizeType>>)
                  {
                    return feeval_ptr->get_symmetric_gradient(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match vector FEEvaluation"));
                    return T {};
                  }
              }
          },
          feeval_vector[(global_variable_index * max_dependency_types) +
                        static_cast<Types::Index>(dependency_type)]);
      }
  }

  /**
   * @brief Return the curl of the specified field.
   *
   * Note that this is dealii::VectorizedArray<number> type for 2D and dealii::Tensor<1,
   * dim, dealii::VectorizedArray<number>> type for 3D.
   *
   * @tparam T the return type. Must be either a `dealii::Tensor<1, 1,
   * SizeType>` or `dealii::Tensor<1, dim, SizeType>`.
   * @param global_variable_index The global index of the variable to access.
   * @param dependency_type The dependency type of the variable to access.
   */
  template <typename T>
  [[nodiscard]] T
  get_vector_curl(Types::Index   global_variable_index,
                  DependencyType dependency_type = DependencyType::Normal) const
  requires(std::is_same_v<T, dealii::Tensor<1, 1, SizeType>> ||
           std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::gradients);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    // For the 1D case the FEEvaluation objects for scalar/vector fields are degenerate.
    // Otherwise, use std::visit
    if constexpr (dim == 1)
      {
        Assert(false, dealii::ExcMessage("Curl is nonsensical for 1D."));
        return T {};
      }
    else
      {
        return std::visit<T>(
          [&](auto &feeval_ptr) -> T
          {
            using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;

            if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
              {
                Assert(false,
                       dealii::ExcMessage("You cannot take the curl of scalar fields."));
                return T {};
              }
            else if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
              {
                if constexpr (std::is_same_v<
                                T,
                                dealii::Tensor<1, (dim == 2 ? 1 : dim), SizeType>>)
                  {
                    return feeval_ptr->get_curl(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Requested type does not match vector FEEvaluation"));
                    return T {};
                  }
              }
          },
          feeval_vector[(global_variable_index * max_dependency_types) +
                        static_cast<Types::Index>(dependency_type)]);
      }
  }

  /**
   * @brief Set the residual value of the specified scalar/vector field.
   */
  template <typename T>
  void
  set_value_term(Types::Index   global_variable_index,
                 const T       &val,
                 DependencyType dependency_type = DependencyType::Normal)
  requires(std::is_same_v<T, SizeType> ||
           std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
  {
    // Some checks that make sure that we have a valid submission and that the
    // feevaluation exists
#ifdef DEBUG
    submission_valid(dependency_type);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    // For the 1D case the FEEvaluation objects for scalar/vector fields are degenerate.
    // Otherwise, use std::visit
    if constexpr (dim == 1)
      {
        // Unwrap the vector value here
        if constexpr (std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
          {
            return feeval_vector[(global_variable_index * max_dependency_types) +
                                 static_cast<Types::Index>(dependency_type)]
              ->submit_value(val[0], q_point);
          }
        return feeval_vector[(global_variable_index * max_dependency_types) +
                             static_cast<Types::Index>(dependency_type)]
          ->submit_value(val, q_point);
      }
    else
      {
        std::visit<void>(
          Overload {[&](std::unique_ptr<ScalarFEEvaluation> &feeval_ptr)
                    {
                      if constexpr (std::is_same_v<T, SizeType>)
                        {
                          feeval_ptr->submit_value(val, q_point);
                        }
                      else
                        {
                          Assert(false,
                                 dealii::ExcMessage(
                                   "Submitted type does not match the one expected for a "
                                   "scalar FEEvaluation object."));
                        }
                    },
                    [&](std::unique_ptr<VectorFEEvaluation> &feeval_ptr)
                    {
                      if constexpr (std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
                        {
                          feeval_ptr->submit_value(val, q_point);
                        }
                      else
                        {
                          Assert(false,
                                 dealii::ExcMessage(
                                   "Submitted type does not match the one expected for a "
                                   "vector FEEvaluation object."));
                        }
                    }},
          feeval_vector[(global_variable_index * max_dependency_types) +
                        static_cast<Types::Index>(dependency_type)]);
      }
  }

  /**
   * @brief Set the residual gradient of the specified scalar/vector field.
   */
  template <typename T>
  void
  set_gradient_term(Types::Index   global_variable_index,
                    const T       &grad,
                    DependencyType dependency_type = DependencyType::Normal)
  requires(std::is_same_v<T, dealii::Tensor<1, dim, SizeType>> ||
           std::is_same_v<T, dealii::Tensor<2, dim, SizeType>>)
  {
    // Some checks that make sure that we have a valid submission and that the
    // feevaluation exists
#ifdef DEBUG
    submission_valid(dependency_type);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    // For the 1D case the FEEvaluation objects for scalar/vector fields are degenerate.
    // Otherwise, use std::visit
    if constexpr (dim == 1)
      {
        // Unwrap the vector value here
        if constexpr (std::is_same_v<T, dealii::Tensor<2, dim, SizeType>>)
          {
            return feeval_vector[(global_variable_index * max_dependency_types) +
                                 static_cast<Types::Index>(dependency_type)]
              ->submit_gradient(grad[0], q_point);
          }
        return feeval_vector[(global_variable_index * max_dependency_types) +
                             static_cast<Types::Index>(dependency_type)]
          ->submit_gradient(grad, q_point);
      }
    else
      {
        std::visit<void>(
          Overload {[&](std::unique_ptr<ScalarFEEvaluation> &feeval_ptr)
                    {
                      if constexpr (std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
                        {
                          feeval_ptr->submit_gradient(grad, q_point);
                        }
                      else
                        {
                          Assert(false,
                                 dealii::ExcMessage(
                                   "Submitted type does not match the one expected for a "
                                   "scalar FEEvaluation object."));
                        }
                    },
                    [&](std::unique_ptr<VectorFEEvaluation> &feeval_ptr)
                    {
                      if constexpr (std::is_same_v<T, dealii::Tensor<2, dim, SizeType>>)
                        {
                          feeval_ptr->submit_gradient(grad, q_point);
                        }
                      else
                        {
                          Assert(false,
                                 dealii::ExcMessage(
                                   "Submitted type does not match the one expected for a "
                                   "vector FEEvaluation object."));
                        }
                    }},
          feeval_vector[(global_variable_index * max_dependency_types) +
                        static_cast<Types::Index>(dependency_type)]);
      }
  }

  /**
   * @brief Apply some operator function for a given cell range and source vector to
   * some destination vector.
   */
  void
  eval_local_operator(const std::function<void(VariableContainer &,
                                               const dealii::Point<dim, SizeType> &,
                                               const SizeType &)> &func,
                      std::vector<VectorType *>                   &dst,
                      const std::vector<VectorType *>             &src,
                      const std::pair<unsigned int, unsigned int> &cell_range);

  /**
   * @brief Apply some operator function for a given cell range and source vector to
   * some destination vector.
   */
  void
  eval_local_operator(const std::function<void(VariableContainer &,
                                               const dealii::Point<dim, SizeType> &,
                                               const SizeType &)> &func,
                      VectorType                                  &dst,
                      const std::vector<VectorType *>             &src,
                      const std::pair<unsigned int, unsigned int> &cell_range);

  /**
   * @brief Apply some operator function for a given cell range and source vector to
   * some destination vector.
   */
  void
  eval_local_operator(const std::function<void(VariableContainer &,
                                               const dealii::Point<dim, SizeType> &,
                                               const SizeType &)> &func,
                      VectorType                                  &dst,
                      const VectorType                            &src,
                      const std::vector<VectorType *>             &src_subset,
                      const std::pair<unsigned int, unsigned int> &cell_range);

  /**
   * @brief TODO (landinjm): Add comments
   */
  void
  eval_local_diagonal(const std::function<void(VariableContainer &,
                                               const dealii::Point<dim, SizeType> &,
                                               const SizeType &)> &func,
                      VectorType                                  &dst,
                      const std::vector<VectorType *>             &src_subset,
                      const std::pair<unsigned int, unsigned int> &cell_range);

private:
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
  using ScalarDiagonal = dealii::AlignedVector<SizeType>;

  /**
   * @brief Typedef for vector diagonal matrix objects.
   */
  using VectorDiagonal = dealii::AlignedVector<dealii::Tensor<1, dim, SizeType>>;

  /**
   * @brief Typedef for the variant diagonal matrix objects.
   */
  using VariantDiagonal =
    std::variant<std::unique_ptr<ScalarDiagonal>, std::unique_ptr<VectorDiagonal>>;

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
  [[nodiscard]] dealii::Point<dim, SizeType>
  get_q_point_location() const;

  /**
   * @brief Initialize, read DOFs, and set evaulation flags for each variable.
   */
  void
  reinit_and_eval(const std::vector<VectorType *> &src, unsigned int cell);

  /**
   * @brief Initialize, read DOFs, and set evaulation flags for each variable.
   */
  void
  reinit_and_eval(const VectorType &src, unsigned int cell);

  /**
   * @brief Initialize the cell for all dependencies of a certain variable index.
   */
  void
  reinit(unsigned int cell, Types::Index global_variable_index);

  /**
   * @brief Read dofs values on the cell for all dependencies of a certain variable index.
   */
  void
  read_dof_values(const std::vector<VectorType *> &src);

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
  integrate_and_distribute(std::vector<VectorType *> &dst);

  /**
   * @brief Integrate the residuals and distribute from local to global.
   */
  void
  integrate_and_distribute(VectorType &dst);

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
  eval_cell_diagonal(FEEvaluationType                            *feeval_ptr,
                     DiagonalType                                *diagonal_ptr,
                     unsigned int                                 cell,
                     Types::Index                                 global_var_index,
                     const std::function<void(VariableContainer &,
                                              const dealii::Point<dim, SizeType> &,
                                              const SizeType &)> &func,
                     VectorType                                  &dst,
                     const std::vector<VectorType *>             &src_subset);

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
