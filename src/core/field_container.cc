// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/exceptions.h>
#include <prismspf/core/field_container.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/utilities/element_volume.h>

#include <prismspf/config.h>

#include "prismspf/core/dependencies.h"
#include "prismspf/core/solution_indexer.h"
#include "prismspf/core/solve_group.h"

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <ranges>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
FieldContainer<dim, degree, number>::FieldContainer(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
  const std::vector<FieldAttributes>       &_field_attributes,
  const SolveGroup                         &_solve_group,
  SolutionIndexer<dim, number>             &_solution_indexer,
  const ElementVolume<dim, degree, number> &_element_volume,
  const std::vector<Types::Index>          &_global_to_local_solution,
  const EquationType                        _equation_type,
  bool                                      use_local_mapping)
  : element_volume_handler(&_element_volume)
  , global_to_local_solution(&_global_to_local_solution)
  , equation_type(equation_type)
{
  const std::vector<FieldAttributes> field_attributes;
  // Initialize the feeval_vector
  feeval_deps_scalar.clear();
  feeval_deps_scalar.resize(field_attributes.size());
  const DependencySet &dependency_map = equation_type == EquationType::RHS
                                          ? solve_group->dependencies_rhs
                                          : solve_group->dependencies_lhs;
  for (const auto &[field_index, dependency] : dependency_map)
    {
      const auto mf_id_pair =
        solution_indexer->get_matrix_free_and_block_index(field_index, relative_level);
      if (field_attributes[field_index].field_type == FieldInfo::TensorRank::Scalar)
        {
          feeval_deps_scalar[field_index] = {dependency, mf_id_pair};
        }
      else if (field_attributes[field_index].field_type == FieldInfo::TensorRank::Vector)
        {
          feeval_deps_vector[field_index] = {dependency, mf_id_pair};
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::eval_local_diagonal(
  const std::function<void(FieldContainer &,
                           const dealii::Point<dim, ScalarValue> &,
                           const ScalarValue &)> &func,
  SolutionVector                                 &dst,
  const std::vector<SolutionVector *>            &src_subset,
  const std::pair<unsigned int, unsigned int>    &cell_range)
{
  Assert(subset_attributes->size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  const auto &global_var_index = subset_attributes->begin()->first;
  const auto &field_type = subset_attributes->begin()->second.field_info.tensor_rank;
  feevaluation_exists(global_var_index, DependencyType::Change);
  auto &feeval_variant = feeval_vector[(global_var_index * max_dependency_types) +
                                       static_cast<Types::Index>(DependencyType::Change)];

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

  if (field_type == FieldInfo::TensorRank::Scalar)
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
  else if (field_type == FieldInfo::TensorRank::Vector)
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
FieldContainer<dim, degree, number>::feevaluation_exists(
  [[maybe_unused]] Types::Index field_index,
  [[maybe_unused]] Types::Index dependency_index) const
{
#ifdef DEBUG
  feevaluation_size_valid(field_index, dependency_index);

  // Check if the feeval_variant is nullptr
  bool is_nullptr = true;
  if constexpr (dim == 1)
    {
      is_nullptr =
        feeval_vector[(field_index * max_dependency_types) + dependency_index] == nullptr;
    }
  else
    {
      is_nullptr = std::visit(
        [](const auto &ptr) -> bool
        {
          return ptr == nullptr;
        },
        feeval_vector[(field_index * max_dependency_types) + dependency_index]);
    }
  Assert(!is_nullptr,
         dealii::ExcMessage("The FEEvaluation object with global index = " +
                            std::to_string(field_index) + " does not exist for type = " +
                            to_string(static_cast<DependencyType>(dependency_index))));
#endif
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::access_valid(
  [[maybe_unused]] Types::Index                             field_index,
  [[maybe_unused]] DependencyType                           dependency_type,
  [[maybe_unused]] dealii::EvaluationFlags::EvaluationFlags flag) const
{
#ifdef DEBUG
  for (const auto &[index, variable] : *subset_attributes)
    {
      if (solve_type == SolveType::NonexplicitLHS)
        {
          Assert(variable.get_eval_flag_set_lhs()[field_index][static_cast<Types::Index>(
                   dependency_type)] != dealii::EvaluationFlags::EvaluationFlags::nothing,
                 DependencyNotFound(field_index, to_string(dependency_type)));
        }
      else
        {
          Assert(variable.get_eval_flag_set_rhs()[field_index][static_cast<Types::Index>(
                   dependency_type)] != dealii::EvaluationFlags::EvaluationFlags::nothing,
                 DependencyNotFound(field_index, to_string(dependency_type)));
        }
    }
#endif
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::submission_valid(
  [[maybe_unused]] DependencyType dependency_type) const
{
#ifdef DEBUG
  Assert(dependency_type == Normal || solve_type == SolveType::NonexplicitLHS,
         dealii::ExcMessage(
           "RHS residuals are only allowed to submit normal gradient terms."));
  Assert(dependency_type == Change || solve_type != SolveType::NonexplicitLHS,
         dealii::ExcMessage(
           "LHS residuals are only allowed to submit change gradient terms."));
#endif
}

template <unsigned int dim, unsigned int degree, typename number>
unsigned int
FieldContainer<dim, degree, number>::get_n_q_points() const
{
  // For dim = 1, the scalar and vector FEEvaluation objects are degenerate.
  if constexpr (dim == 1)
    {
      auto iterator = std::ranges::find_if(feeval_vector,
                                           [](const auto &ptr)
                                           {
                                             return ptr != nullptr;
                                           });
      Assert(iterator != feeval_vector.end(),
             dealii::ExcMessage("All FEEvaluation objects were nullptr."));
      return (*iterator)->n_q_points;
    }
  else
    {
      // Create a filtered view that ignores nullptrs
      auto filtered_view = feeval_vector | std::views::filter(
                                             [](const auto &v)
                                             {
                                               return !std::visit(
                                                 [](const auto &ptr)
                                                 {
                                                   return ptr == nullptr;
                                                 },
                                                 v);
                                             });
      // Grab the first iterator and return the number of quad points
      auto iterator = filtered_view.begin();
      Assert(iterator != filtered_view.end(),
             dealii::ExcMessage("All FEEvaluation variants were nullptr."));
      return std::visit(
        [](const auto &ptr)
        {
          return ptr->n_q_points;
        },
        *iterator);
    }
  Assert(false,
         dealii::ExcMessage("When trying to access the number of quadrature points, all "
                            "FEEvaluation object containers were empty."));
  return 0;
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Point<dim, typename FieldContainer<dim, degree, number>::ScalarValue>
FieldContainer<dim, degree, number>::get_q_point_location() const
{
  // For dim = 1, the scalar and vector FEEvaluation objects are degenerate.
  if constexpr (dim == 1)
    {
      auto iterator = std::ranges::find_if(feeval_vector,
                                           [](const auto &ptr)
                                           {
                                             return ptr != nullptr;
                                           });
      Assert(iterator != feeval_vector.end(),
             dealii::ExcMessage("All FEEvaluation objects were nullptr."));
      return (*iterator)->quadrature_point(q_point);
    }
  else
    {
      // Create a filtered view that ignores nullptrs
      auto filtered_view = feeval_vector | std::views::filter(
                                             [](const auto &v)
                                             {
                                               return !std::visit(
                                                 [](const auto &ptr)
                                                 {
                                                   return ptr == nullptr;
                                                 },
                                                 v);
                                             });
      // Grab the first iterator and return the number of quad points
      auto iterator = filtered_view.begin();
      Assert(iterator != filtered_view.end(),
             dealii::ExcMessage("All FEEvaluation variants were nullptr."));
      return std::visit(
        [this](const auto &ptr)
        {
          return ptr->quadrature_point(q_point);
        },
        *iterator);
    }
  Assert(false,
         dealii::ExcMessage("When trying to access the quadrature point location, all "
                            "FEEvaluation object containers were empty."));
  return dealii::Point<dim, ScalarValue>();
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::reinit(unsigned int cell)
{
  const DependencySet &dependency_map; // rhs or lhs
  for (const auto &[field_index, dependency] : dependency_map)
    {
      feeval_deps_scalar[field_index].reinit(cell);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::eval()
{
  const DependencySet &dependency_map; // rhs or lhs
  for (const auto &[field_index, dependency] : dependency_map)
    {
      feeval_deps_scalar[field_index].eval();
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::reinit_and_eval(unsigned int cell)
{
  const DependencySet &dependency_map; // rhs or lhs
  for (const auto &[field_index, dependency] : dependency_map)
    {
      feeval_deps_scalar[field_index].reinit_and_eval(cell);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::integrate(Types::Index global_variable_index)
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
        feeval_vector[(global_variable_index * max_dependency_types) +
                      static_cast<Types::Index>(DependencyType::Change)];

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
                            "Unexpected type in feeval_vector variant");
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
FieldContainer<dim, degree, number>::integrate_and_distribute(SolutionVector &dst)
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
    auto &feeval_variant = feeval_vector[(residual_index * max_dependency_types) +
                                         static_cast<Types::Index>(dependency_type)];

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
                          "Unexpected type in feeval_vector variant");
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
template <typename FEEvaluationType, typename DiagonalType>
void
FieldContainer<dim, degree, number>::eval_cell_diagonal(
  FEEvaluationType                               *feeval_ptr,
  DiagonalType                                   *diagonal_ptr,
  unsigned int                                    cell,
  Types::Index                                    global_variable_index,
  const std::function<void(FieldContainer &,
                           const dealii::Point<dim, ScalarValue> &,
                           const ScalarValue &)> &func,
  SolutionVector                                 &dst,
  const std::vector<SolutionVector *>            &src_subset)
{
  using DiagonalValueType = typename DiagonalType::value_type;

  // Grab the element volume
  const ScalarValue element_volume = element_volume_handler->get_element_volume(cell);

  // Reinit the cell for all the dependencies
  reinit(cell, global_variable_index);

  for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
    {
      int dof_index = i;
      // Submit an identity matrix for the change term
      for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
        {
          if constexpr (std::is_same_v<DiagonalValueType, ScalarValue> || dim == 1)
            {
              feeval_ptr->submit_dof_value(ScalarValue(), j);
            }
          else
            {
              feeval_ptr->submit_dof_value(DiagonalValueType(), j);
            }
        }

      // Set the i-th value to 1.0
      if constexpr (std::is_same_v<DiagonalValueType, ScalarValue> || dim == 1)
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

      // Read plain dof values for non change src
      read_dof_values(src_subset);

      // Evaluate the dependencies based on the flags
      eval(global_variable_index);

      // Evaluate at each quadrature point
      for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
        {
          q_point = quad;
          func(*this, get_q_point_location(), element_volume);
        }

      // Integrate the diagonal
      integrate(global_variable_index);
      (*diagonal_ptr)[i] = feeval_ptr->get_dof_value(i);
    }

  // Submit calculated diagonal values and distribute
  for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
    {
      if constexpr (std::is_same_v<DiagonalValueType, ScalarValue> || dim != 1)
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

#include "core/variable_container.inst"

PRISMS_PF_END_NAMESPACE
