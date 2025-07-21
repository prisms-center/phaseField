// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/exceptions.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/config.h>

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
VariableContainer<dim, degree, number>::VariableContainer(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
  const std::map<unsigned int, VariableAttributes> &_subset_attributes,
  const std::vector<std::vector<Types::Index>>     &_global_to_local_solution,
  const SolveType                                  &_solve_type,
  bool                                              use_local_mapping)
  : subset_attributes(&_subset_attributes)
  , global_to_local_solution(&_global_to_local_solution)
  , solve_type(_solve_type)
{
  // Initialize the feeval_map
  const Types::Index max_fields =
    subset_attributes->begin()->second.get_dependency_set_rhs().size();
  const Types::Index max_dependencies =
    subset_attributes->begin()->second.get_dependency_set_rhs().begin()->size();
  feeval_map.resize(max_fields);
  for (auto &dependency_feeval_map : feeval_map)
    {
      dependency_feeval_map.resize(max_dependencies);
    }

  auto construct_map = [&](const std::vector<std::vector<FieldType>> &dependency_set)
  {
    Types::Index dependency_index = 0;
    for (const auto &inner_dependency_set : dependency_set)
      {
        Types::Index dependency_type = 0;
        for (const auto &field_type : inner_dependency_set)
          {
            // Skip if the field type is invalid
            if (field_type == Numbers::invalid_field_type)
              {
                dependency_type++;
                continue;
              }

            // Ensure the indices are within bounds
            Assert(dependency_index < feeval_map.size(),
                   dealii::ExcMessage(
                     "Dependency index " + std::to_string(dependency_index) +
                     " exceeds feeval_map size " + std::to_string(feeval_map.size())));
            Assert(dependency_type < feeval_map[dependency_index].size(),
                   dealii::ExcMessage(
                     "Dependency type index " +
                     to_string(static_cast<DependencyType>(dependency_type)) +
                     " exceeds feeval_map[dependency_index] size " +
                     std::to_string(feeval_map[dependency_index].size())));

            if (field_type == FieldType::Scalar)
              {
                if (!use_local_mapping)
                  {
                    feeval_map[dependency_index][dependency_type] =
                      std::make_unique<ScalarFEEvaluation>(data, dependency_index);
                  }
                else
                  {
                    // TODO (landinjm): Find a better way to represent this. Maybe it
                    // would be better just to pass the desired mapping of global
                    // variables to matrix free indices. For most cases, they are one and
                    // the same, but for multigrid they are different as not all fields
                    // have matrixfree data associated for the multigrid levels.

                    Assert(global_to_local_solution->size() > dependency_index,
                           dealii::ExcMessage(
                             "Global to local solution size " +
                             std::to_string(global_to_local_solution->size()) +
                             " is greater than dependency index " +
                             std::to_string(dependency_index)));
                    Assert(global_to_local_solution->at(dependency_index).size() >
                             dependency_type,
                           dealii::ExcMessage(
                             "Global to local solution size at dependency index " +
                             std::to_string(dependency_index) + " " +
                             std::to_string(
                               global_to_local_solution->at(dependency_index).size()) +
                             " is greater than dependency type " +
                             std::to_string(dependency_type)));
                    Assert(global_to_local_solution->at(dependency_index)
                               .at(dependency_type) != Numbers::invalid_index,
                           dealii::ExcMessage(
                             "Global to local solution at dependency index " +
                             std::to_string(dependency_index) + " " +
                             std::to_string(dependency_type) + " is invalid"));

                    feeval_map[dependency_index][dependency_type] =
                      std::make_unique<ScalarFEEvaluation>(data,
                                                           global_to_local_solution
                                                             ->at(dependency_index)
                                                             .at(dependency_type));
                  }
              }
            else
              {
                if (!use_local_mapping)
                  {
                    feeval_map[dependency_index][dependency_type] =
                      std::make_unique<VectorFEEvaluation>(data, dependency_index);
                  }
                else
                  {
                    // TODO (landinjm): Find a better way to represent this. Maybe it
                    // would be better just to pass the desired mapping of global
                    // variables to matrix free indices. For most cases, they are one and
                    // the same, but for multigrid they are different as not all fields
                    // have matrixfree data associated for the multigrid levels.

                    Assert(global_to_local_solution->size() > dependency_index,
                           dealii::ExcMessage(
                             "Global to local solution size " +
                             std::to_string(global_to_local_solution->size()) +
                             " is greater than dependency index " +
                             std::to_string(dependency_index)));
                    Assert(global_to_local_solution->at(dependency_index).size() >
                             dependency_type,
                           dealii::ExcMessage(
                             "Global to local solution size at dependency index " +
                             std::to_string(dependency_index) + " " +
                             std::to_string(
                               global_to_local_solution->at(dependency_index).size()) +
                             " is greater than dependency type " +
                             std::to_string(dependency_type)));
                    Assert(global_to_local_solution->at(dependency_index)
                               .at(dependency_type) != Numbers::invalid_index,
                           dealii::ExcMessage(
                             "Global to local solution at dependency index " +
                             std::to_string(dependency_index) + " " +
                             std::to_string(dependency_type) + " is invalid"));

                    feeval_map[dependency_index][dependency_type] =
                      std::make_unique<VectorFEEvaluation>(data,
                                                           global_to_local_solution
                                                             ->at(dependency_index)
                                                             .at(dependency_type));
                  }
              }

            dependency_type++;
          }

        dependency_index++;
      }
  };

  Assert(subset_attributes->size() == 1 ||
           (solve_type == SolveType::ExplicitRHS || solve_type == SolveType::Postprocess),
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  if (solve_type == SolveType::NonexplicitLHS)
    {
      construct_map(subset_attributes->begin()->second.get_dependency_set_lhs());
      return;
    }
  construct_map(subset_attributes->begin()->second.get_dependency_set_rhs());
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::eval_local_operator(
  const std::function<void(VariableContainer &, const dealii::Point<dim, SizeType> &)>
                                              &func,
  std::vector<VectorType *>                   &dst,
  const std::vector<VectorType *>             &src,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      reinit_and_eval(src, cell);

      // Evaluate at each quadrature point
      for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
        {
          q_point = quad;
          func(*this, get_q_point_location());
        }

      // Integrate and add to global vector dst
      integrate_and_distribute(dst);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::eval_local_operator(
  const std::function<void(VariableContainer &, const dealii::Point<dim, SizeType> &)>
                                              &func,
  VectorType                                  &dst,
  const std::vector<VectorType *>             &src,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      reinit_and_eval(src, cell);

      // Evaluate at each quadrature point
      for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
        {
          q_point = quad;
          func(*this, get_q_point_location());
        }

      // Integrate and add to global vector dst
      integrate_and_distribute(dst);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::eval_local_operator(
  const std::function<void(VariableContainer &, const dealii::Point<dim, SizeType> &)>
                                              &func,
  VectorType                                  &dst,
  const VectorType                            &src,
  const std::vector<VectorType *>             &src_subset,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      reinit_and_eval(src, cell);
      reinit_and_eval(src_subset, cell);

      // Evaluate at each quadrature point
      for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
        {
          q_point = quad;
          func(*this, get_q_point_location());
        }

      // Integrate and add to global vector dst
      integrate_and_distribute(dst);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
template <typename FEEvaluationType>
FEEvaluationType *
VariableContainer<dim, degree, number>::extract_feeval_ptr(
  VariantFEEvaluation &variant) const
{
  FEEvaluationType *return_ptr = nullptr;
  std::visit(
    [&return_ptr](auto &ptr)
    {
      using T = std::decay_t<decltype(ptr)>;
      if constexpr (std::is_same_v<T, std::unique_ptr<FEEvaluationType>>)
        {
          return_ptr = ptr.get();
        }
      else
        {
          Assert(false, dealii::ExcNotInitialized());
        }
    },
    variant);
  return return_ptr;
}

template <unsigned int dim, unsigned int degree, typename number>
template <typename DiagonalType>
DiagonalType *
VariableContainer<dim, degree, number>::extract_diagonal_ptr(
  VariantDiagonal &variant) const
{
  DiagonalType *return_ptr = nullptr;
  std::visit(
    [&return_ptr](auto &ptr)
    {
      using T = std::decay_t<decltype(ptr)>;
      if constexpr (std::is_same_v<T, std::unique_ptr<DiagonalType>>)
        {
          return_ptr = ptr.get();
        }
      else
        {
          Assert(false, dealii::ExcNotInitialized());
        }
    },
    variant);
  return return_ptr;
}

template <unsigned int dim, unsigned int degree, typename number>
template <typename FEEvaluationType, typename DiagonalType>
void
VariableContainer<dim, degree, number>::eval_cell_diagonal(
  FEEvaluationType *feeval_ptr,
  DiagonalType     *diagonal_ptr,
  unsigned int      cell,
  unsigned int      global_var_index,
  const std::function<void(VariableContainer &, const dealii::Point<dim, SizeType> &)>
                                  &func,
  VectorType                      &dst,
  const std::vector<VectorType *> &src_subset)
{
  using DiagonalValueType = typename DiagonalType::value_type;

  // Helper function to submit the identity matrix
  auto submit_identity = [&](auto &feeval_ptr, unsigned int dof_index)
  {
    for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
      {
        if constexpr (std::is_same_v<DiagonalValueType, SizeType> || dim == 1)
          {
            feeval_ptr->submit_dof_value(SizeType(), j);
          }
        else
          {
            feeval_ptr->submit_dof_value(DiagonalValueType(), j);
          }
      }

    // Set the i-th value to 1.0
    if constexpr (std::is_same_v<DiagonalValueType, SizeType> || dim == 1)
      {
        feeval_ptr->submit_dof_value(dealii::make_vectorized_array<number>(1.0),
                                     dof_index);
      }
    else
      {
        DiagonalValueType one;
        for (unsigned int dimension = 0; dimension < dim; ++dimension)
          {
            one[dimension] = dealii::make_vectorized_array<number>(1.0);
          }
        feeval_ptr->submit_dof_value(one, dof_index);
      }
  };
  // Reinit the cell for all the dependencies
  reinit(cell, global_var_index);

  for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
    {
      // Submit an identity matrix for the change term
      submit_identity(feeval_ptr, i);

      // Read plain dof values for non change src
      read_dof_values(src_subset);

      // Evaluate the dependencies based on the flags
      eval(global_var_index);

      // Evaluate at each quadrature point
      for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
        {
          q_point = quad;
          func(*this, get_q_point_location());
        }

      // Integrate the diagonal
      integrate(global_var_index);
      (*diagonal_ptr)[i] = feeval_ptr->get_dof_value(i);
    }

  // Submit calculated diagonal values and distribute
  for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
    {
      if constexpr (std::is_same_v<DiagonalValueType, SizeType> || dim != 1)
        {
          feeval_ptr->submit_dof_value((*diagonal_ptr)[i], i);
        }
      else
        {
          feeval_ptr->submit_dof_value((*diagonal_ptr)[i][0], i);
        }
    }
  feeval_ptr->distribute_local_to_global(dst);
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::eval_local_diagonal(
  const std::function<void(VariableContainer &, const dealii::Point<dim, SizeType> &)>
                                              &func,
  VectorType                                  &dst,
  const std::vector<VectorType *>             &src_subset,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  Assert(subset_attributes->size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  const auto &global_var_index = subset_attributes->begin()->first;
  const auto &field_type       = subset_attributes->begin()->second.get_field_type();
  feevaluation_exists(global_var_index, DependencyType::Change);
  auto &feeval_variant =
    feeval_map[global_var_index][static_cast<Types::Index>(DependencyType::Change)];

  auto process_feeval = [&](auto &feeval_ptr, auto &diag_ptr)
  {
    using FEEvaluationType = std::decay_t<decltype(*feeval_ptr)>;
    using DiagonalType     = std::decay_t<decltype(*diag_ptr)>;

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        eval_cell_diagonal<FEEvaluationType, DiagonalType>(feeval_ptr,
                                                           diag_ptr,
                                                           cell,
                                                           global_var_index,
                                                           func,
                                                           dst,
                                                           src_subset);
      }
  };

  if (field_type == FieldType::Scalar)
    {
      ScalarFEEvaluation *scalar_feeval_ptr = nullptr;

      if constexpr (dim == 1)
        {
          scalar_feeval_ptr = feeval_variant.get();
        }
      else
        {
          scalar_feeval_ptr = extract_feeval_ptr<ScalarFEEvaluation>(feeval_variant);
        }

      n_dofs_per_cell       = scalar_feeval_ptr->dofs_per_cell;
      diagonal              = std::make_unique<ScalarDiagonal>(n_dofs_per_cell);
      auto *scalar_diag_ptr = extract_diagonal_ptr<ScalarDiagonal>(diagonal);
      process_feeval(scalar_feeval_ptr, scalar_diag_ptr);
    }
  else if (field_type == FieldType::Vector)
    {
      VectorFEEvaluation *vector_feeval_ptr = nullptr;

      if constexpr (dim == 1)
        {
          vector_feeval_ptr = feeval_variant.get();
        }
      else
        {
          vector_feeval_ptr = extract_feeval_ptr<VectorFEEvaluation>(feeval_variant);
        }

      n_dofs_per_cell       = vector_feeval_ptr->dofs_per_component;
      diagonal              = std::make_unique<VectorDiagonal>(n_dofs_per_cell);
      auto *vector_diag_ptr = extract_diagonal_ptr<VectorDiagonal>(diagonal);
      process_feeval(vector_feeval_ptr, vector_diag_ptr);
    }
  else
    {
      Assert(false, UnreachableCode());
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::feevaluation_exists(
  [[maybe_unused]] const unsigned int   &dependency_index,
  [[maybe_unused]] const DependencyType &dependency_type) const
{
#ifdef DEBUG
  Assert(feeval_map.size() > dependency_index,
         dealii::ExcMessage("The FEEvaluation object does not exist for global index = " +
                            std::to_string(dependency_index)));
  Assert(feeval_map.at(dependency_index).size() >
           static_cast<Types::Index>(dependency_type),
         dealii::ExcMessage("The FEEvaluation object with global index = " +
                            std::to_string(dependency_index) +
                            " does not exist for type = " + to_string(dependency_type)));
  // Check if the feeval_variant is nullptr
  bool is_nullptr = true;

  if constexpr (dim == 1)
    {
      is_nullptr =
        feeval_map[dependency_index][static_cast<Types::Index>(dependency_type)] ==
        nullptr;
    }
  else
    {
      is_nullptr = std::visit(
        [](const auto &ptr) -> bool
        {
          return ptr == nullptr;
        },
        feeval_map[dependency_index][static_cast<Types::Index>(dependency_type)]);
    }

  Assert(!is_nullptr,
         dealii::ExcMessage("The FEEvaluation object with global index = " +
                            std::to_string(dependency_index) +
                            " does not exist for type = " + to_string(dependency_type)));
#endif
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::access_valid(
  [[maybe_unused]] const unsigned int                             &dependency_index,
  [[maybe_unused]] const DependencyType                           &dependency_type,
  [[maybe_unused]] const dealii::EvaluationFlags::EvaluationFlags &flag) const
{
  for ([[maybe_unused]] const auto &[index, variable] : *subset_attributes)
    {
      if (solve_type == SolveType::NonexplicitLHS)
        {
          Assert(
            variable.get_eval_flag_set_lhs()[dependency_index][static_cast<Types::Index>(
              dependency_type)] != dealii::EvaluationFlags::EvaluationFlags::nothing,
            DependencyNotFound(dependency_index, to_string(dependency_type)));
        }
      else
        {
          Assert(
            variable.get_eval_flag_set_rhs()[dependency_index][static_cast<Types::Index>(
              dependency_type)] != dealii::EvaluationFlags::EvaluationFlags::nothing,
            DependencyNotFound(dependency_index, to_string(dependency_type)));
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::submission_valid(
  [[maybe_unused]] const DependencyType &dependency_type) const
{
  Assert(dependency_type == Normal || solve_type == SolveType::NonexplicitLHS,
         dealii::ExcMessage(
           "RHS residuals are only allowed to submit normal gradient terms."));
  Assert(dependency_type == Change || solve_type != SolveType::NonexplicitLHS,
         dealii::ExcMessage(
           "LHS residuals are only allowed to submit change gradient terms."));
}

template <unsigned int dim, unsigned int degree, typename number>
unsigned int
VariableContainer<dim, degree, number>::get_n_q_points() const
{
  for (const auto &dependency_feeval_map : feeval_map)
    {
      for (const auto &feeval_variant : dependency_feeval_map)
        {
          if constexpr (dim == 1)
            {
              if (feeval_variant == nullptr)
                {
                  continue;
                }
              return feeval_variant->n_q_points;
            }
          else
            {
              bool is_nullptr = std::visit(
                [](const auto &ptr) -> bool
                {
                  return ptr == nullptr;
                },
                feeval_variant);

              if (is_nullptr)
                {
                  continue;
                }

              return std::visit(
                [&](const auto &feeval_ptr) -> unsigned int
                {
                  Assert(feeval_ptr != nullptr, dealii::ExcNotInitialized());
                  return feeval_ptr->n_q_points;
                },
                feeval_variant);
            }
        }
    }

  Assert(false,
         dealii::ExcMessage("When trying to access the number of quadrature points, all "
                            "FEEvaluation object containers were empty."));
  return 0;
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Point<dim, typename VariableContainer<dim, degree, number>::SizeType>
VariableContainer<dim, degree, number>::get_q_point_location() const
{
  for (const auto &dependency_feeval_map : feeval_map)
    {
      for (const auto &feeval_variant : dependency_feeval_map)
        {
          if constexpr (dim == 1)
            {
              if (feeval_variant == nullptr)
                {
                  continue;
                }
              return feeval_variant->quadrature_point(q_point);
            }
          else
            {
              bool is_nullptr = std::visit(
                [](const auto &ptr) -> bool
                {
                  return ptr == nullptr;
                },
                feeval_variant);

              if (is_nullptr)
                {
                  continue;
                }

              return std::visit(
                [&](const auto &feeval_ptr) -> dealii::Point<dim, SizeType>
                {
                  Assert(feeval_ptr != nullptr, dealii::ExcNotInitialized());
                  return feeval_ptr->quadrature_point(q_point);
                },
                feeval_variant);
            }
        }
    }

  Assert(false,
         dealii::ExcMessage("When trying to access the quadrature point location, all "
                            "FEEvaluation object containers were empty."));

  return dealii::Point<dim, SizeType>();
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::reinit_and_eval(
  const std::vector<VectorType *> &src,
  unsigned int                     cell)
{
  // Reinit and eval values for the given dependency set. Note the dependency set includes
  // the variable we're evaluating, which may or may not be an actually dependency. For
  // this reason, I selectively read dofs and evaluate the flags.
  auto reinit_and_eval_map =
    [&](const std::vector<std::vector<dealii::EvaluationFlags::EvaluationFlags>>
                                                  &eval_flag_set,
        const std::vector<std::vector<FieldType>> &dependency_set)
  {
    Types::Index dependency_index = 0;
    for (const auto &inner_dependency_set : dependency_set)
      {
        Types::Index dependency_type = 0;
        for (const auto &field_type : inner_dependency_set)
          {
            if (static_cast<DependencyType>(dependency_type) == DependencyType::Change ||
                field_type == Numbers::invalid_field_type)
              {
                dependency_type++;
                continue;
              }

            Assert(global_to_local_solution->size() > dependency_index,
                   dealii::ExcMessage("Global to local solution size " +
                                      std::to_string(global_to_local_solution->size()) +
                                      " is greater than dependency index " +
                                      std::to_string(dependency_index)));
            Assert(
              global_to_local_solution->at(dependency_index).size() > dependency_type,
              dealii::ExcMessage(
                "Global to local solution size at dependency index " +
                std::to_string(dependency_index) + " " +
                std::to_string(global_to_local_solution->at(dependency_index).size()) +
                " is greater than dependency type " + std::to_string(dependency_type)));
            Assert(global_to_local_solution->at(dependency_index).at(dependency_type) !=
                     Numbers::invalid_index,
                   dealii::ExcMessage("Global to local solution at dependency index " +
                                      std::to_string(dependency_index) + " " +
                                      std::to_string(dependency_type) + " is invalid"));
            feevaluation_exists(dependency_index,
                                static_cast<DependencyType>(dependency_type));
            auto &feeval_variant = feeval_map[dependency_index][dependency_type];

            auto process_feeval = [&](auto &feeval_ptr)
            {
              feeval_ptr->reinit(cell);
              if (eval_flag_set[dependency_index][dependency_type] ==
                  dealii::EvaluationFlags::EvaluationFlags::nothing)
                {
                  return;
                }
              const Types::Index &local_index =
                global_to_local_solution->at(dependency_index).at(dependency_type);
              Assert(src.size() > local_index,
                     dealii::ExcMessage(
                       "The provided src vector's size is below the given local "
                       "index = " +
                       std::to_string(local_index) + " for global index = " +
                       std::to_string(dependency_index) + "  and type = " +
                       to_string(static_cast<DependencyType>(dependency_type))));
              feeval_ptr->read_dof_values_plain(*(src.at(local_index)));
              feeval_ptr->evaluate(eval_flag_set[dependency_index][dependency_type]);
            };

            if constexpr (dim == 1)
              {
                process_feeval(feeval_variant);
              }
            else
              {
                std::visit(
                  [&](auto &ptr)
                  {
                    using T = std::decay_t<decltype(ptr)>;
                    static_assert(
                      std::is_same_v<T, std::unique_ptr<ScalarFEEvaluation>> ||
                        std::is_same_v<T, std::unique_ptr<VectorFEEvaluation>>,
                      "Unexpected type in feeval_map variant");
                    process_feeval(ptr);
                  },
                  feeval_variant);
              }

            dependency_type++;
          }

        dependency_index++;
      }
  };

  if (solve_type == SolveType::ExplicitRHS || solve_type == SolveType::Postprocess)
    {
      const auto &attrs = subset_attributes->begin()->second;
      reinit_and_eval_map(attrs.get_eval_flag_set_rhs(), attrs.get_dependency_set_rhs());
      return;
    }
  if (src.empty())
    {
      return;
    }

  Assert(subset_attributes->size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));
  const auto &attrs = subset_attributes->begin()->second;
  if (solve_type == SolveType::NonexplicitLHS)
    {
      reinit_and_eval_map(attrs.get_eval_flag_set_lhs(), attrs.get_dependency_set_lhs());
    }
  else
    {
      reinit_and_eval_map(attrs.get_eval_flag_set_rhs(), attrs.get_dependency_set_rhs());
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::reinit_and_eval(const VectorType &src,
                                                        unsigned int      cell)
{
  // Reinit and eval values for the given dependency set. Note the dependency set includes
  // the variable we're evaluating, which may or may not be an actually dependency. For
  // this reason, I selectively read dofs and evaluate the flags.
  auto reinit_and_eval_map =
    [&](const std::vector<std::vector<dealii::EvaluationFlags::EvaluationFlags>>
                                                  &eval_flag_set,
        const std::vector<std::vector<FieldType>> &dependency_set)
  {
    Types::Index dependency_index = 0;
    for (const auto &inner_dependency_set : dependency_set)
      {
        Types::Index dependency_type = 0;
        for (const auto &field_type : inner_dependency_set)
          {
            // TODO (landinjm): This can be drastically simplified because all we're doing
            // in reinit-ing and eval-ing the change solution.
            if (static_cast<DependencyType>(dependency_type) != DependencyType::Change ||
                field_type == Numbers::invalid_field_type)
              {
                dependency_type++;
                continue;
              }

            feevaluation_exists(dependency_index,
                                static_cast<DependencyType>(dependency_type));
            auto &feeval_variant = feeval_map[dependency_index][dependency_type];

            auto process_feeval = [&](auto &feeval_ptr)
            {
              feeval_ptr->reinit(cell);
              feeval_ptr->read_dof_values_plain(src);
              feeval_ptr->evaluate(eval_flag_set[dependency_index][dependency_type]);
            };

            if constexpr (dim == 1)
              {
                process_feeval(feeval_variant);
              }
            else
              {
                std::visit(
                  [&](auto &ptr)
                  {
                    using T = std::decay_t<decltype(ptr)>;
                    static_assert(
                      std::is_same_v<T, std::unique_ptr<ScalarFEEvaluation>> ||
                        std::is_same_v<T, std::unique_ptr<VectorFEEvaluation>>,
                      "Unexpected type in feeval_map variant");
                    process_feeval(ptr);
                  },
                  feeval_variant);
              }

            dependency_type++;
          }

        dependency_index++;
      }
  };

  Assert(subset_attributes->size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));
  const auto &attrs = subset_attributes->begin()->second;
  if (solve_type == SolveType::NonexplicitLHS)
    {
      reinit_and_eval_map(attrs.get_eval_flag_set_lhs(), attrs.get_dependency_set_lhs());
      return;
    }
  Assert(false,
         dealii::ExcMessage(
           "reinit_and_eval(src) should only be called for LHS evaluations"));
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::reinit(unsigned int        cell,
                                               const unsigned int &global_variable_index)
{
  auto reinit_map =
    [&](const std::vector<std::vector<dealii::EvaluationFlags::EvaluationFlags>>
          &eval_flag_set)
  {
    Types::Index dependency_index = 0;
    for (const auto &dependency_set : eval_flag_set)
      {
        Types::Index dependency_type = 0;
        for (const auto &value : dependency_set)
          {
            if (value == dealii::EvaluationFlags::EvaluationFlags::nothing)
              {
                dependency_type++;
                continue;
              }

            feevaluation_exists(dependency_index,
                                static_cast<DependencyType>(dependency_type));
            auto &feeval_variant =
              feeval_map[dependency_index][static_cast<Types::Index>(dependency_type)];

            auto process_feeval = [&](auto &feeval_ptr)
            {
              feeval_ptr->reinit(cell);
            };

            if constexpr (dim == 1)
              {
                process_feeval(feeval_variant);
              }
            else
              {
                std::visit(
                  [&](auto &ptr)
                  {
                    using T = std::decay_t<decltype(ptr)>;
                    static_assert(
                      std::is_same_v<T, std::unique_ptr<ScalarFEEvaluation>> ||
                        std::is_same_v<T, std::unique_ptr<VectorFEEvaluation>>,
                      "Unexpected type in feeval_map variant");
                    process_feeval(ptr);
                  },
                  feeval_variant);
              }
            dependency_type++;
          }

        dependency_index++;
      }
  };

  Assert(subset_attributes->contains(global_variable_index),
         dealii::ExcMessage(
           "The subset attribute entry does not exists for global index = " +
           std::to_string(global_variable_index)));
  const auto &attrs = subset_attributes->at(global_variable_index);
  if (solve_type == SolveType::NonexplicitLHS)
    {
      reinit_map(attrs.get_eval_flag_set_lhs());
    }
  else
    {
      reinit_map(attrs.get_eval_flag_set_rhs());
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::read_dof_values(
  const std::vector<VectorType *> &src)
{
  auto reinit_and_eval_map =
    [&](const std::vector<std::vector<dealii::EvaluationFlags::EvaluationFlags>>
                                                  &eval_flag_set,
        const std::vector<std::vector<FieldType>> &dependency_set)
  {
    Types::Index dependency_index = 0;
    for (const auto &inner_dependency_set : dependency_set)
      {
        Types::Index dependency_type = 0;
        for (const auto &field_type : inner_dependency_set)
          {
            if (static_cast<DependencyType>(dependency_type) == DependencyType::Change ||
                field_type == Numbers::invalid_field_type)
              {
                dependency_type++;
                continue;
              }

            Assert(global_to_local_solution->size() > dependency_index,
                   dealii::ExcMessage("Global to local solution size " +
                                      std::to_string(global_to_local_solution->size()) +
                                      " is greater than dependency index " +
                                      std::to_string(dependency_index)));
            Assert(
              global_to_local_solution->at(dependency_index).size() > dependency_type,
              dealii::ExcMessage(
                "Global to local solution size at dependency index " +
                std::to_string(dependency_index) + " " +
                std::to_string(global_to_local_solution->at(dependency_index).size()) +
                " is greater than dependency type " + std::to_string(dependency_type)));
            Assert(global_to_local_solution->at(dependency_index).at(dependency_type) !=
                     Numbers::invalid_index,
                   dealii::ExcMessage("Global to local solution at dependency index " +
                                      std::to_string(dependency_index) + " " +
                                      std::to_string(dependency_type) + " is invalid"));
            feevaluation_exists(dependency_index,
                                static_cast<DependencyType>(dependency_type));
            auto &feeval_variant = feeval_map[dependency_index][dependency_type];

            auto process_feeval = [&](auto &feeval_ptr)
            {
              if (eval_flag_set[dependency_index][dependency_type] !=
                  dealii::EvaluationFlags::EvaluationFlags::nothing)
                {
                  const Types::Index &local_index =
                    global_to_local_solution->at(dependency_index).at(dependency_type);
                  Assert(src.size() > local_index,
                         dealii::ExcMessage(
                           "The provided src vector's size is below the given local "
                           "index = " +
                           std::to_string(local_index) + " for global index = " +
                           std::to_string(dependency_index) + "  and type = " +
                           to_string(static_cast<DependencyType>(dependency_type))));
                  feeval_ptr->read_dof_values_plain(*(src.at(local_index)));
                }
            };

            if constexpr (dim == 1)
              {
                process_feeval(feeval_variant);
              }
            else
              {
                std::visit(
                  [&](auto &ptr)
                  {
                    using T = std::decay_t<decltype(ptr)>;
                    static_assert(
                      std::is_same_v<T, std::unique_ptr<ScalarFEEvaluation>> ||
                        std::is_same_v<T, std::unique_ptr<VectorFEEvaluation>>,
                      "Unexpected type in feeval_map variant");
                    process_feeval(ptr);
                  },
                  feeval_variant);
              }

            dependency_type++;
          }

        dependency_index++;
      }
  };

  if (solve_type != SolveType::NonexplicitLHS)
    {
      Assert(false, UnreachableCode());
      return;
    }
  if (src.empty())
    {
      return;
    }

  Assert(subset_attributes->size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));
  const auto &attrs = subset_attributes->begin()->second;
  reinit_and_eval_map(attrs.get_eval_flag_set_lhs(), attrs.get_dependency_set_lhs());
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::eval(const unsigned int &global_variable_index)
{
  auto eval_map =
    [&](const std::vector<std::vector<dealii::EvaluationFlags::EvaluationFlags>>
          &eval_flag_set)
  {
    Types::Index index = 0;
    for (const auto &dependency_set : eval_flag_set)
      {
        Types::Index dep_index = 0;
        for (const auto &value : dependency_set)
          {
            if (value == dealii::EvaluationFlags::EvaluationFlags::nothing)
              {
                dep_index++;
                continue;
              }

            feevaluation_exists(index, static_cast<DependencyType>(dep_index));
            auto &feeval_variant = feeval_map[index][dep_index];

            auto process_feeval = [&](auto &feeval_ptr)
            {
              feeval_ptr->evaluate(value);
            };

            if constexpr (dim == 1)
              {
                process_feeval(feeval_variant);
              }
            else
              {
                std::visit(
                  [&](auto &ptr)
                  {
                    using T = std::decay_t<decltype(ptr)>;
                    static_assert(
                      std::is_same_v<T, std::unique_ptr<ScalarFEEvaluation>> ||
                        std::is_same_v<T, std::unique_ptr<VectorFEEvaluation>>,
                      "Unexpected type in feeval_map variant");
                    process_feeval(ptr);
                  },
                  feeval_variant);
              }

            dep_index++;
          }

        index++;
      }
  };

  Assert(subset_attributes->contains(global_variable_index),
         dealii::ExcMessage(
           "The subset attribute entry does not exists for global index = " +
           std::to_string(global_variable_index)));
  const auto &attrs = subset_attributes->at(global_variable_index);
  if (solve_type == SolveType::NonexplicitLHS)
    {
      eval_map(attrs.get_eval_flag_set_lhs());
    }
  else
    {
      eval_map(attrs.get_eval_flag_set_rhs());
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::integrate(
  const unsigned int &global_variable_index)
{
  Assert(subset_attributes->contains(global_variable_index),
         dealii::ExcMessage(
           "The subset attribute entry does not exists for global index = " +
           std::to_string(global_variable_index)));
  const auto &variable = subset_attributes->at(global_variable_index);

  if (solve_type == SolveType::NonexplicitLHS)
    {
      feevaluation_exists(global_variable_index, DependencyType::Change);
      auto &feeval_variant =
        feeval_map[global_variable_index]
                  [static_cast<Types::Index>(DependencyType::Change)];

      auto process_feeval = [&](auto &feeval_ptr)
      {
        feeval_ptr->integrate(variable.get_eval_flags_residual_lhs());
      };

      if constexpr (dim == 1)
        {
          process_feeval(feeval_variant);
        }
      else
        {
          std::visit(
            [&](auto &ptr)
            {
              using T = std::decay_t<decltype(ptr)>;
              static_assert(std::is_same_v<T, std::unique_ptr<ScalarFEEvaluation>> ||
                              std::is_same_v<T, std::unique_ptr<VectorFEEvaluation>>,
                            "Unexpected type in feeval_map variant");
              process_feeval(ptr);
            },
            feeval_variant);
        }
      return;
    }

  AssertThrow(false,
              dealii::ExcMessage(
                "Integrate called for a solve type that is not NonexplicitLHS."));
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::integrate_and_distribute(
  std::vector<VectorType *> &dst)
{
  auto integrate_and_distribute_map =
    [&](const dealii::EvaluationFlags::EvaluationFlags &residual_flag_set,
        const DependencyType                           &dependency_type,
        const unsigned int                             &residual_index)
  {
    Assert(global_to_local_solution->size() > residual_index,
           dealii::ExcMessage("Global to local solution size " +
                              std::to_string(global_to_local_solution->size()) +
                              " is greater than dependency index " +
                              std::to_string(residual_index)));
    Assert(global_to_local_solution->at(residual_index).size() > dependency_type,
           dealii::ExcMessage(
             "Global to local solution size at dependency index " +
             std::to_string(residual_index) + " " +
             std::to_string(global_to_local_solution->at(residual_index).size()) +
             " is greater than dependency type " + std::to_string(dependency_type)));
    Assert(global_to_local_solution->at(residual_index).at(dependency_type) !=
             Numbers::invalid_index,
           dealii::ExcMessage("Global to local solution at dependency index " +
                              std::to_string(residual_index) + " " +
                              std::to_string(dependency_type) + " is invalid"));
    const Types::Index &local_index =
      global_to_local_solution->at(residual_index).at(dependency_type);
    Assert(dst.size() > local_index,
           dealii::ExcMessage(
             "The provided dst vector's size is below the given local index = " +
             std::to_string(local_index) +
             " for global index = " + std::to_string(residual_index) +
             "  and type = " + to_string(dependency_type)));
    Assert(subset_attributes->contains(residual_index),
           dealii::ExcMessage(
             "The subset attribute entry does not exists for global index = " +
             std::to_string(residual_index)));
    feevaluation_exists(residual_index, dependency_type);
    auto &feeval_variant =
      feeval_map[residual_index][static_cast<Types::Index>(dependency_type)];

    auto process_feeval = [&](auto &feeval_ptr)
    {
      feeval_ptr->integrate_scatter(residual_flag_set, *(dst.at(local_index)));
    };

    if constexpr (dim == 1)
      {
        process_feeval(feeval_variant);
      }
    else
      {
        std::visit(
          [&](auto &ptr)
          {
            using T = std::decay_t<decltype(ptr)>;
            static_assert(std::is_same_v<T, std::unique_ptr<ScalarFEEvaluation>> ||
                            std::is_same_v<T, std::unique_ptr<VectorFEEvaluation>>,
                          "Unexpected type in feeval_map variant");
            process_feeval(ptr);
          },
          feeval_variant);
      }
  };

  for (const auto &[index, variable] : *subset_attributes)
    {
      if (solve_type == SolveType::NonexplicitLHS)
        {
          integrate_and_distribute_map(variable.get_eval_flags_residual_lhs(),
                                       DependencyType::Change,
                                       index);
        }
      else
        {
          integrate_and_distribute_map(variable.get_eval_flags_residual_rhs(),
                                       DependencyType::Normal,
                                       index);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::integrate_and_distribute(VectorType &dst)
{
  auto integrate_and_distribute_map =
    [&](const dealii::EvaluationFlags::EvaluationFlags &residual_flag_set,
        const DependencyType                           &dependency_type,
        const unsigned int                             &residual_index)
  {
    Assert(subset_attributes->contains(residual_index),
           dealii::ExcMessage(
             "The subset attribute entry does not exists for global index = " +
             std::to_string(residual_index)));
    feevaluation_exists(residual_index, dependency_type);
    auto &feeval_variant =
      feeval_map[residual_index][static_cast<Types::Index>(dependency_type)];

    auto process_feeval = [&](auto &feeval_ptr)
    {
      feeval_ptr->integrate_scatter(residual_flag_set, dst);
    };

    if constexpr (dim == 1)
      {
        process_feeval(feeval_variant);
      }
    else
      {
        std::visit(
          [&](auto &ptr)
          {
            using T = std::decay_t<decltype(ptr)>;
            static_assert(std::is_same_v<T, std::unique_ptr<ScalarFEEvaluation>> ||
                            std::is_same_v<T, std::unique_ptr<VectorFEEvaluation>>,
                          "Unexpected type in feeval_map variant");
            process_feeval(ptr);
          },
          feeval_variant);
      }
  };

  Assert(subset_attributes->size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  if (solve_type == SolveType::NonexplicitLHS)
    {
      integrate_and_distribute_map(
        subset_attributes->begin()->second.get_eval_flags_residual_lhs(),
        DependencyType::Change,
        subset_attributes->begin()->first);
    }
  else
    {
      integrate_and_distribute_map(
        subset_attributes->begin()->second.get_eval_flags_residual_rhs(),
        DependencyType::Normal,
        subset_attributes->begin()->first);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
typename VariableContainer<dim, degree, number>::SizeType
VariableContainer<dim, degree, number>::get_vector_divergence(
  unsigned int   global_variable_index,
  DependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      return feeval_map[global_variable_index][static_cast<Types::Index>(dependency_type)]
        ->get_divergence(q_point);
    }
  else
    {
      return std::visit(
        [&](const auto &feeval_ptr) -> SizeType
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
            {
              return feeval_ptr->get_divergence(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage(
                       "Expected VectorFEEvaluation but got ScalarFEEvaluation."));
              return SizeType();
            }
        },
        feeval_map[global_variable_index][static_cast<Types::Index>(dependency_type)]);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Tensor<2, dim, typename VariableContainer<dim, degree, number>::SizeType>
VariableContainer<dim, degree, number>::get_vector_symmetric_gradient(
  unsigned int   global_variable_index,
  DependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      return feeval_map[global_variable_index][static_cast<Types::Index>(dependency_type)]
        ->get_symmetric_gradient(q_point);
    }
  else
    {
      return std::visit(
        [&](const auto &feeval_ptr) -> dealii::Tensor<2, dim, SizeType>
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
            {
              return feeval_ptr->get_symmetric_gradient(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage(
                       "Expected VectorFEEvaluation but got ScalarFEEvaluation."));
              return dealii::Tensor<2, dim, SizeType>();
            }
        },
        feeval_map[global_variable_index][static_cast<Types::Index>(dependency_type)]);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Tensor<1,
               (dim == 2 ? 1 : dim),
               typename VariableContainer<dim, degree, number>::SizeType>
VariableContainer<dim, degree, number>::get_vector_curl(
  unsigned int   global_variable_index,
  DependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      Assert(false, dealii::ExcMessage("Curl is nonsensical for 1D."));
      return dealii::Tensor<1, (dim == 2 ? 1 : dim), SizeType> {};
    }
  else
    {
      // Return the value directly for dim > 1
      return std::visit(
        [&](const auto &feeval_ptr) -> dealii::Tensor<1, (dim == 2 ? 1 : dim), SizeType>
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
            {
              return feeval_ptr->get_curl(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage(
                       "Expected VectorFEEvaluation but got ScalarFEEvaluation."));
              return dealii::Tensor<1, (dim == 2 ? 1 : dim), SizeType>();
            }
        },
        feeval_map[global_variable_index][static_cast<Types::Index>(dependency_type)]);
    }
}

INSTANTIATE_TRI_TEMPLATE(VariableContainer)

PRISMS_PF_END_NAMESPACE
