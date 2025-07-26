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
    const std::map<unsigned int, VariableAttributes> &_subset_attributes,
    const std::vector<std::vector<Types::Index>>     &_global_to_local_solution,
    const SolveType                                  &_solve_type,
    bool                                              use_local_mapping = false);

  /**
   * @brief Return the value of the specified field.
   */
  template <FieldType T>
  [[nodiscard]] std::
    conditional_t<T == FieldType::Scalar, SizeType, dealii::Tensor<1, dim, SizeType>>
    get_value(Types::Index   global_variable_index,
              DependencyType dependency_type = DependencyType::Normal) const
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::values);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    if constexpr (T == FieldType::Scalar)
      {
        if constexpr (dim == 1)
          {
            return feeval_map[global_variable_index]
                             [static_cast<Types::Index>(dependency_type)]
                               ->get_value(q_point);
          }
        else
          {
            return std::visit(
              [&](const auto &feeval_ptr) -> SizeType
              {
                using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
                if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
                  {
                    return feeval_ptr->get_value(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Expected ScalarFEEvaluation but got VectorFEEvaluation."));
                    return SizeType();
                  }
              },
              feeval_map[global_variable_index]
                        [static_cast<Types::Index>(dependency_type)]);
          }
      }
    else
      {
        if constexpr (dim == 1)
          {
            dealii::Tensor<1, dim, SizeType> wrapper;
            wrapper[0] = feeval_map[global_variable_index]
                                   [static_cast<Types::Index>(dependency_type)]
                                     ->get_value(q_point);
            return wrapper;
          }
        else
          {
            return std::visit(
              [&](const auto &feeval_ptr) -> dealii::Tensor<1, dim, SizeType>
              {
                using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
                if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
                  {
                    return feeval_ptr->get_value(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Expected VectorFEEvaluation but got ScalarFEEvaluation."));
                    return dealii::Tensor<1, dim, SizeType>();
                  }
              },
              feeval_map[global_variable_index]
                        [static_cast<Types::Index>(dependency_type)]);
          }
      }
  }

  /**
   * @brief Return the gradient of the specified field.
   */
  template <FieldType T>
  [[nodiscard]] std::conditional_t<T == FieldType::Scalar,
                                   dealii::Tensor<1, dim, SizeType>,
                                   dealii::Tensor<2, dim, SizeType>>
  get_gradient(Types::Index   global_variable_index,
               DependencyType dependency_type = DependencyType::Normal) const
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::gradients);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    if constexpr (T == FieldType::Scalar)
      {
        if constexpr (dim == 1)
          {
            return feeval_map[global_variable_index]
                             [static_cast<Types::Index>(dependency_type)]
                               ->get_gradient(q_point);
          }
        else
          {
            return std::visit(
              [&](const auto &feeval_ptr) -> dealii::Tensor<1, dim, SizeType>
              {
                using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
                if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
                  {
                    return feeval_ptr->get_gradient(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Expected ScalarFEEvaluation but got VectorFEEvaluation."));
                    return dealii::Tensor<1, dim, SizeType>();
                  }
              },
              feeval_map[global_variable_index]
                        [static_cast<Types::Index>(dependency_type)]);
          }
      }
    else
      {
        if constexpr (dim == 1)
          {
            dealii::Tensor<2, dim, SizeType> wrapper;
            wrapper[0] = feeval_map[global_variable_index]
                                   [static_cast<Types::Index>(dependency_type)]
                                     ->get_gradient(q_point);
            return wrapper;
          }
        else
          {
            return std::visit(
              [&](const auto &feeval_ptr) -> dealii::Tensor<2, dim, SizeType>
              {
                using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
                if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
                  {
                    return feeval_ptr->get_gradient(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Expected VectorFEEvaluation but got ScalarFEEvaluation."));
                    return dealii::Tensor<2, dim, SizeType>();
                  }
              },
              feeval_map[global_variable_index]
                        [static_cast<Types::Index>(dependency_type)]);
          }
      }
  }

  /**
   * @brief Return the hessian of the specified field.
   */
  template <FieldType T>
  [[nodiscard]] std::conditional_t<T == FieldType::Scalar,
                                   dealii::Tensor<2, dim, SizeType>,
                                   dealii::Tensor<3, dim, SizeType>>
  get_hessian(Types::Index   global_variable_index,
              DependencyType dependency_type = DependencyType::Normal) const
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::hessians);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    if constexpr (T == FieldType::Scalar)
      {
        if constexpr (dim == 1)
          {
            return feeval_map[global_variable_index]
                             [static_cast<Types::Index>(dependency_type)]
                               ->get_hessian(q_point);
          }
        else
          {
            return std::visit(
              [&](const auto &feeval_ptr) -> dealii::Tensor<2, dim, SizeType>
              {
                using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
                if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
                  {
                    return feeval_ptr->get_hessian(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Expected ScalarFEEvaluation but got VectorFEEvaluation."));
                    return dealii::Tensor<2, dim, SizeType>();
                  }
              },
              feeval_map[global_variable_index]
                        [static_cast<Types::Index>(dependency_type)]);
          }
      }
    else
      {
        if constexpr (dim == 1)
          {
            dealii::Tensor<3, dim, SizeType> wrapper;
            wrapper[0] = feeval_map[global_variable_index]
                                   [static_cast<Types::Index>(dependency_type)]
                                     ->get_hessian(q_point);
            return wrapper;
          }
        else
          {
            return std::visit(
              [&](const auto &feeval_ptr) -> dealii::Tensor<3, dim, SizeType>
              {
                using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
                if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
                  {
                    return feeval_ptr->get_hessian(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Expected VectorFEEvaluation but got ScalarFEEvaluation."));
                    return dealii::Tensor<3, dim, SizeType>();
                  }
              },
              feeval_map[global_variable_index]
                        [static_cast<Types::Index>(dependency_type)]);
          }
      }
  }

  /**
   * @brief Return the diagonal of the hessian of the specified field.
   */
  template <FieldType T>
  [[nodiscard]] std::conditional_t<T == FieldType::Scalar,
                                   dealii::Tensor<1, dim, SizeType>,
                                   dealii::Tensor<2, dim, SizeType>>
  get_hessian_diagonal(Types::Index   global_variable_index,
                       DependencyType dependency_type = DependencyType::Normal) const
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::hessians);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    if constexpr (T == FieldType::Scalar)
      {
        if constexpr (dim == 1)
          {
            return feeval_map[global_variable_index]
                             [static_cast<Types::Index>(dependency_type)]
                               ->get_hessian_diagonal(q_point);
          }
        else
          {
            return std::visit(
              [&](const auto &feeval_ptr) -> dealii::Tensor<1, dim, SizeType>
              {
                using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
                if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
                  {
                    return feeval_ptr->get_hessian_diagonal(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Expected ScalarFEEvaluation but got VectorFEEvaluation."));
                    return dealii::Tensor<1, dim, SizeType>();
                  }
              },
              feeval_map[global_variable_index]
                        [static_cast<Types::Index>(dependency_type)]);
          }
      }
    else
      {
        if constexpr (dim == 1)
          {
            dealii::Tensor<2, dim, SizeType> wrapper;
            wrapper[0] = feeval_map[global_variable_index]
                                   [static_cast<Types::Index>(dependency_type)]
                                     ->get_hessian_diagonal(q_point);
            return wrapper;
          }
        else
          {
            return std::visit(
              [&](const auto &feeval_ptr) -> dealii::Tensor<2, dim, SizeType>
              {
                using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
                if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
                  {
                    return feeval_ptr->get_hessian_diagonal(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Expected VectorFEEvaluation but got ScalarFEEvaluation."));
                    return dealii::Tensor<2, dim, SizeType>();
                  }
              },
              feeval_map[global_variable_index]
                        [static_cast<Types::Index>(dependency_type)]);
          }
      }
  }

  /**
   * @brief Return the laplacian of the specified field.
   */
  template <FieldType T>
  [[nodiscard]] std::
    conditional_t<T == FieldType::Scalar, SizeType, dealii::Tensor<1, dim, SizeType>>
    get_laplacian(Types::Index   global_variable_index,
                  DependencyType dependency_type = DependencyType::Normal) const
  {
    // Some checks that make sure that we have evaluation flags and that the feevaluation
    // exists
#ifdef DEBUG
    access_valid(global_variable_index,
                 dependency_type,
                 dealii::EvaluationFlags::EvaluationFlags::hessians);
    feevaluation_exists(global_variable_index, dependency_type);
#endif

    if constexpr (T == FieldType::Scalar)
      {
        if constexpr (dim == 1)
          {
            return feeval_map[global_variable_index]
                             [static_cast<Types::Index>(dependency_type)]
                               ->get_laplacian(q_point);
          }
        else
          {
            return std::visit(
              [&](const auto &feeval_ptr) -> SizeType
              {
                using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
                if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
                  {
                    return feeval_ptr->get_laplacian(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Expected ScalarFEEvaluation but got VectorFEEvaluation."));
                    return SizeType();
                  }
              },
              feeval_map[global_variable_index]
                        [static_cast<Types::Index>(dependency_type)]);
          }
      }
    else
      {
        if constexpr (dim == 1)
          {
            dealii::Tensor<1, dim, SizeType> wrapper;
            wrapper[0] = feeval_map[global_variable_index]
                                   [static_cast<Types::Index>(dependency_type)]
                                     ->get_laplacian(q_point);
            return wrapper;
          }
        else
          {
            return std::visit(
              [&](const auto &feeval_ptr) -> dealii::Tensor<1, dim, SizeType>
              {
                using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
                if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
                  {
                    return feeval_ptr->get_laplacian(q_point);
                  }
                else
                  {
                    Assert(false,
                           dealii::ExcMessage(
                             "Expected VectorFEEvaluation but got ScalarFEEvaluation."));
                    return dealii::Tensor<1, dim, SizeType>();
                  }
              },
              feeval_map[global_variable_index]
                        [static_cast<Types::Index>(dependency_type)]);
          }
      }
  }

  /**
   * @brief Return the divergence of the specified vector field.
   *
   * TODO (landinjm): Not sure if we should bother template this to make it look like the
   * other function. It will just throw an error in the scalar case.
   */
  [[nodiscard]] SizeType
  get_vector_divergence(unsigned int   global_variable_index,
                        DependencyType dependency_type = DependencyType::Normal) const;

  /**
   * @brief Return the symmetric gradient of the specified vector field.
   */
  [[nodiscard]] dealii::Tensor<2, dim, SizeType>
  get_vector_symmetric_gradient(
    unsigned int   global_variable_index,
    DependencyType dependency_type = DependencyType::Normal) const;

  /**
   * @brief Return the curl of the specified vector field. Note that this is
   * dealii::VectorizedArray<number> type for 2D and dealii::Tensor<1, dim,
   * dealii::VectorizedArray<number>> type for 3D.
   */
  [[nodiscard]] dealii::Tensor<1, (dim == 2 ? 1 : dim), SizeType>
  get_vector_curl(unsigned int   global_variable_index,
                  DependencyType dependency_type = DependencyType::Normal) const;

  /**
   * @brief Set the residual value of the specified scalar/vector field.
   */
  template <typename T>
  void
  set_value_term(const unsigned int   &global_variable_index,
                 const T              &val,
                 const DependencyType &dependency_type = DependencyType::Normal)
  requires(std::is_same_v<T, SizeType> ||
           std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>)
  {
    // Some checks that make sure that we have a valid submission and that the
    // feevaluation exists
    // TODO (Landinjm): Check that the feevaluation is the right type
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
            return feeval_map[global_variable_index]
                             [static_cast<Types::Index>(dependency_type)]
                               ->submit_value(val[0], q_point);
          }
        return feeval_map[global_variable_index]
                         [static_cast<Types::Index>(dependency_type)]
                           ->submit_value(val, q_point);
      }
    else
      {
        std::visit(
          [&](auto &feeval_ptr)
          {
            using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
            if constexpr (FEEvalType::n_components == 1)
              {
                static_assert(std::is_same_v<T, SizeType>,
                              "Expected SizeType for scalar field.");
                feeval_ptr->submit_value(val, q_point);
              }
            else if constexpr (FEEvalType::n_components == dim)
              {
                static_assert(
                  std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>,
                  "Expected dealii::Tensor<1, dim, SizeType> for vector field.");
                feeval_ptr->submit_value(val, q_point);
              }
            else
              {
                static_assert(FEEvalType::n_components == 1 ||
                                FEEvalType::n_components == dim,
                              "Unexpected number of components");
              }
          },
          feeval_map[global_variable_index][static_cast<Types::Index>(dependency_type)]);
      }
  }

  /**
   * @brief Set the residual gradient of the specified scalar/vector field.
   */
  template <typename T>
  void
  set_gradient_term(const unsigned int   &global_variable_index,
                    const T              &grad,
                    const DependencyType &dependency_type = DependencyType::Normal)
  requires(std::is_same_v<T, dealii::Tensor<1, dim, SizeType>> ||
           std::is_same_v<T, dealii::Tensor<2, dim, SizeType>>)
  {
    // Some checks that make sure that we have a valid submission and that the
    // feevaluation exists
    // TODO (Landinjm): Check that the feevaluation is the right type
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
            return feeval_map[global_variable_index]
                             [static_cast<Types::Index>(dependency_type)]
                               ->submit_gradient(grad[0], q_point);
          }
        return feeval_map[global_variable_index]
                         [static_cast<Types::Index>(dependency_type)]
                           ->submit_gradient(grad, q_point);
      }
    else
      {
        std::visit(
          [&](auto &feeval_ptr)
          {
            using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
            if constexpr (FEEvalType::n_components == 1)
              {
                static_assert(
                  std::is_same_v<T, dealii::Tensor<1, dim, SizeType>>,
                  "Expected dealii::Tensor<1, dim, SizeType> for scalar field.");
                feeval_ptr->submit_gradient(grad, q_point);
              }
            else if constexpr (FEEvalType::n_components == dim)
              {
                static_assert(
                  std::is_same_v<T, dealii::Tensor<2, dim, SizeType>>,
                  "Expected dealii::Tensor<2, dim, SizeType> for vector field.");
                feeval_ptr->submit_gradient(grad, q_point);
              }
            else
              {
                static_assert(FEEvalType::n_components == 1 ||
                                FEEvalType::n_components == dim,
                              "Unexpected number of components");
              }
          },
          feeval_map[global_variable_index][static_cast<Types::Index>(dependency_type)]);
      }
  }

  /**
   * @brief Apply some operator function for a given cell range and source vector to
   * some destination vector.
   */
  void
  eval_local_operator(
    const std::function<void(VariableContainer &, const dealii::Point<dim, SizeType> &)>
                                                &func,
    std::vector<VectorType *>                   &dst,
    const std::vector<VectorType *>             &src,
    const std::pair<unsigned int, unsigned int> &cell_range);

  /**
   * @brief Apply some operator function for a given cell range and source vector to
   * some destination vector.
   */
  void
  eval_local_operator(
    const std::function<void(VariableContainer &, const dealii::Point<dim, SizeType> &)>
                                                &func,
    VectorType                                  &dst,
    const std::vector<VectorType *>             &src,
    const std::pair<unsigned int, unsigned int> &cell_range);

  /**
   * @brief Apply some operator function for a given cell range and source vector to
   * some destination vector.
   */
  void
  eval_local_operator(
    const std::function<void(VariableContainer &, const dealii::Point<dim, SizeType> &)>
                                                &func,
    VectorType                                  &dst,
    const VectorType                            &src,
    const std::vector<VectorType *>             &src_subset,
    const std::pair<unsigned int, unsigned int> &cell_range);

  /**
   * @brief TODO (landinjm): Add comments
   */
  void
  eval_local_diagonal(
    const std::function<void(VariableContainer &, const dealii::Point<dim, SizeType> &)>
                                                &func,
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
   * @brief Check whether the map entry for the  FEEvaluation exists.
   */
  void
  feevaluation_exists(const unsigned int   &dependency_index,
                      const DependencyType &dependency_type) const;

  /**
   * @brief Check that a variable value/gradient/hessians was marked as needed and thus
   * properly initialized.
   */
  void
  access_valid(const unsigned int                             &dependency_index,
               const DependencyType                           &dependency_type,
               const dealii::EvaluationFlags::EvaluationFlags &flag) const;

  /**
   * @brief Check that a value is valid for submission.
   */
  void
  submission_valid(const DependencyType &dependency_type) const;

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
  reinit(unsigned int cell, const unsigned int &global_variable_index);

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
  eval(const unsigned int &global_variable_index);

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
  integrate(const unsigned int &global_variable_index);

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
  eval_cell_diagonal(
    FEEvaluationType *feeval_ptr,
    DiagonalType     *diagonal_ptr,
    unsigned int      cell,
    unsigned int      global_var_index,
    const std::function<void(VariableContainer &, const dealii::Point<dim, SizeType> &)>
                                    &func,
    VectorType                      &dst,
    const std::vector<VectorType *> &src_subset);

  /**
   * @brief Map of FEEvaluation objects for each active variable. The first mapping is
   * for the global variable, the second is for the DependencyType, and the value is
   * a variant that can hold either a scalar or vector FEEvaluation.
   */
  std::vector<std::vector<VariantFEEvaluation>> feeval_map;

  /**
   * @brief The attribute list of the relevant subset of variables.
   */
  const std::map<unsigned int, VariableAttributes> *subset_attributes;

  /**
   * @brief Mapping from global solution vectors to the local ones
   */
  const std::vector<std::vector<Types::Index>> *global_to_local_solution;

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
