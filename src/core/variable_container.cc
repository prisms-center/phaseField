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

#include <prismspf/utilities/element_volume.h>

#include <prismspf/config.h>

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
VariableContainer<dim, degree, number>::VariableContainer(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
  const std::map<Types::Index, VariableAttributes> &_subset_attributes,
  const ElementVolume<dim, degree, number>         &_element_volume,
  const std::vector<Types::Index>                  &_global_to_local_solution,
  const SolveType                                  &_solve_type,
  bool                                              use_local_mapping)
  : subset_attributes(&_subset_attributes)
  , element_volume_handler(&_element_volume)
  , global_to_local_solution(&_global_to_local_solution)
  , solve_type(_solve_type)
{
  // Grab some data from the VariableAttributes
  max_fields           = subset_attributes->begin()->second.get_max_fields();
  max_dependency_types = subset_attributes->begin()->second.get_max_dependency_types();

  // Initialize the feeval_vector
  feeval_vector.resize(max_fields * max_dependency_types);

  auto construct_map =
    [&](const std::vector<std::vector<FieldInfo::TensorRank>> &dependency_set)
  {
    Types::Index dependency_index = 0;
    for (const auto &inner_dependency_set : dependency_set)
      {
        Types::Index dependency_type = 0;
        for (const auto &field_type : inner_dependency_set)
          {
            // Skip if the field type is invalid
            if (field_type == FieldInfo::TensorRank::Undefined)
              {
                dependency_type++;
                continue;
              }

            feevaluation_size_valid(dependency_index, dependency_type);

            if (field_type == FieldInfo::TensorRank::Scalar)
              {
                if (!use_local_mapping)
                  {
                    feeval_vector[(dependency_index * max_dependency_types) +
                                  dependency_type] =
                      std::make_unique<ScalarFEEvaluation>(data, dependency_index);
                  }
                else
                  {
                    // TODO (landinjm): Find a better way to represent this. Maybe it
                    // would be better just to pass the desired mapping of global
                    // variables to matrix free indices. For most cases, they are one and
                    // the same, but for multigrid they are different as not all fields
                    // have matrixfree data associated for the multigrid levels.
                    const Types::Index local_index =
                      get_local_solution_index(dependency_index, dependency_type);
                    feeval_vector[(dependency_index * max_dependency_types) +
                                  dependency_type] =
                      std::make_unique<ScalarFEEvaluation>(data, local_index);
                  }
              }
            else
              {
                if (!use_local_mapping)
                  {
                    feeval_vector[(dependency_index * max_dependency_types) +
                                  dependency_type] =
                      std::make_unique<VectorFEEvaluation>(data, dependency_index);
                  }
                else
                  {
                    // TODO (landinjm): Find a better way to represent this. Maybe it
                    // would be better just to pass the desired mapping of global
                    // variables to matrix free indices. For most cases, they are one and
                    // the same, but for multigrid they are different as not all fields
                    // have matrixfree data associated for the multigrid levels.
                    const Types::Index local_index =
                      get_local_solution_index(dependency_index, dependency_type);
                    feeval_vector[(dependency_index * max_dependency_types) +
                                  dependency_type] =
                      std::make_unique<VectorFEEvaluation>(data, local_index);
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
  const std::function<void(VariableContainer &,
                           const dealii::Point<dim, SizeType> &,
                           const SizeType &)> &func,
  std::vector<VectorType *>                   &dst,
  const std::vector<VectorType *>             &src,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Grab the element volume
      const SizeType element_volume = element_volume_handler->get_element_volume(cell);

      // Initialize, read DOFs, and set evaulation flags for each variable
      reinit_and_eval(src, cell);

      // Evaluate at each quadrature point
      for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
        {
          q_point = quad;
          func(*this, get_q_point_location(), element_volume);
        }

      // Integrate and add to global vector dst
      integrate_and_distribute(dst);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::eval_local_operator(
  const std::function<void(VariableContainer &,
                           const dealii::Point<dim, SizeType> &,
                           const SizeType &)> &func,
  VectorType                                  &dst,
  const std::vector<VectorType *>             &src,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Grab the element volume
      const SizeType element_volume = element_volume_handler->get_element_volume(cell);

      // Initialize, read DOFs, and set evaulation flags for each variable
      reinit_and_eval(src, cell);

      // Evaluate at each quadrature point
      for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
        {
          q_point = quad;
          func(*this, get_q_point_location(), element_volume);
        }

      // Integrate and add to global vector dst
      integrate_and_distribute(dst);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::eval_local_operator(
  const std::function<void(VariableContainer &,
                           const dealii::Point<dim, SizeType> &,
                           const SizeType &)> &func,
  VectorType                                  &dst,
  const VectorType                            &src,
  const std::vector<VectorType *>             &src_subset,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Grab the element volume
      const SizeType element_volume = element_volume_handler->get_element_volume(cell);

      // Initialize, read DOFs, and set evaulation flags for each variable
      reinit_and_eval(src, cell);
      reinit_and_eval(src_subset, cell);

      // Evaluate at each quadrature point
      for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
        {
          q_point = quad;
          func(*this, get_q_point_location(), element_volume);
        }

      // Integrate and add to global vector dst
      integrate_and_distribute(dst);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::eval_local_diagonal(
  const std::function<void(VariableContainer &,
                           const dealii::Point<dim, SizeType> &,
                           const SizeType &)> &func,
  VectorType                                  &dst,
  const std::vector<VectorType *>             &src_subset,
  const std::pair<unsigned int, unsigned int> &cell_range)
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
VariableContainer<dim, degree, number>::feevaluation_size_valid(
  [[maybe_unused]] Types::Index field_index,
  [[maybe_unused]] Types::Index dependency_index) const
{
#ifdef DEBUG
  Assert(feeval_vector.size() > field_index * max_dependency_types + dependency_index,
         dealii::ExcMessage(
           "The FEEvaluation object with global index = " + std::to_string(field_index) +
           " and type = " + to_string(static_cast<DependencyType>(dependency_index)) +
           " does not exist"));
#endif
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::feevaluation_exists(
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
VariableContainer<dim, degree, number>::global_to_local_solution_exists(
  [[maybe_unused]] Types::Index field_index,
  [[maybe_unused]] Types::Index dependency_index) const
{
#ifdef DEBUG
  Assert(global_to_local_solution->size() >
           field_index * max_dependency_types + dependency_index,
         dealii::ExcMessage(
           "The global to local solution vector size " +
           std::to_string(global_to_local_solution->size()) +
           " is greater than index = " + std::to_string(field_index) +
           " and type = " + to_string(static_cast<DependencyType>(dependency_index))));
  Assert(global_to_local_solution->at((field_index * max_dependency_types) +
                                      dependency_index) != Numbers::invalid_index,
         dealii::ExcMessage("Global to local solution at dependency index " +
                            std::to_string(field_index) + " " +
                            std::to_string(dependency_index) + " is invalid"));

#endif
}

template <unsigned int dim, unsigned int degree, typename number>
Types::Index
VariableContainer<dim, degree, number>::get_local_solution_index(
  Types::Index field_index,
  Types::Index dependency_index) const
{
  // Check that the local index exists and is valid first
#ifdef DEBUG
  global_to_local_solution_exists(field_index, dependency_index);
#endif
  return (
    *global_to_local_solution)[(field_index * max_dependency_types) + dependency_index];
}

template <unsigned int dim, unsigned int degree, typename number>
Types::Index
VariableContainer<dim, degree, number>::get_local_solution_index(
  Types::Index   field_index,
  DependencyType dependency_type) const
{
  return get_local_solution_index(field_index,
                                  static_cast<Types::Index>(dependency_type));
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::access_valid(
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
VariableContainer<dim, degree, number>::submission_valid(
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
VariableContainer<dim, degree, number>::get_n_q_points() const
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
dealii::Point<dim, typename VariableContainer<dim, degree, number>::SizeType>
VariableContainer<dim, degree, number>::get_q_point_location() const
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
        const std::vector<std::vector<FieldInfo::TensorRank>> &dependency_set)
  {
    Types::Index dependency_index = 0;
    for (const auto &inner_dependency_set : dependency_set)
      {
        Types::Index dependency_type = 0;
        for (const auto &field_type : inner_dependency_set)
          {
            if (static_cast<DependencyType>(dependency_type) == DependencyType::Change ||
                field_type == FieldInfo::TensorRank::Undefined)
              {
                dependency_type++;
                continue;
              }

            const Types::Index local_index =
              get_local_solution_index(dependency_index, dependency_type);

            feevaluation_exists(dependency_index,
                                static_cast<DependencyType>(dependency_type));
            auto &feeval_variant =
              feeval_vector[(dependency_index * max_dependency_types) + dependency_type];

            auto process_feeval = [&](auto &feeval_ptr)
            {
              feeval_ptr->reinit(cell);
              if (eval_flag_set[dependency_index][dependency_type] ==
                  dealii::EvaluationFlags::EvaluationFlags::nothing)
                {
                  return;
                }
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
                      "Unexpected type in feeval_vector variant");
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
        const std::vector<std::vector<FieldInfo::TensorRank>> &dependency_set)
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
                field_type == FieldInfo::TensorRank::Undefined)
              {
                dependency_type++;
                continue;
              }

            feevaluation_exists(dependency_index,
                                static_cast<DependencyType>(dependency_type));
            auto &feeval_variant =
              feeval_vector[(dependency_index * max_dependency_types) + dependency_type];

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
                      "Unexpected type in feeval_vector variant");
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
VariableContainer<dim, degree, number>::reinit(unsigned int cell,
                                               Types::Index global_variable_index)
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
              feeval_vector[(dependency_index * max_dependency_types) + dependency_type];

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
                      "Unexpected type in feeval_vector variant");
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
        const std::vector<std::vector<FieldInfo::TensorRank>> &dependency_set)
  {
    Types::Index dependency_index = 0;
    for (const auto &inner_dependency_set : dependency_set)
      {
        Types::Index dependency_type = 0;
        for (const auto &field_type : inner_dependency_set)
          {
            if (static_cast<DependencyType>(dependency_type) == DependencyType::Change ||
                field_type == FieldInfo::TensorRank::Undefined)
              {
                dependency_type++;
                continue;
              }

            const Types::Index local_index =
              get_local_solution_index(dependency_index, dependency_type);

            feevaluation_exists(dependency_index,
                                static_cast<DependencyType>(dependency_type));
            auto &feeval_variant =
              feeval_vector[(dependency_index * max_dependency_types) + dependency_type];

            auto process_feeval = [&](auto &feeval_ptr)
            {
              if (eval_flag_set[dependency_index][dependency_type] !=
                  dealii::EvaluationFlags::EvaluationFlags::nothing)
                {
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
                      "Unexpected type in feeval_vector variant");
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
VariableContainer<dim, degree, number>::eval(Types::Index global_variable_index)
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
            auto &feeval_variant =
              feeval_vector[(index * max_dependency_types) + dep_index];

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
                      "Unexpected type in feeval_vector variant");
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
VariableContainer<dim, degree, number>::integrate(Types::Index global_variable_index)
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
VariableContainer<dim, degree, number>::integrate_and_distribute(
  std::vector<VectorType *> &dst)
{
  auto integrate_and_distribute_map =
    [&](const dealii::EvaluationFlags::EvaluationFlags &residual_flag_set,
        const DependencyType                           &dependency_type,
        const unsigned int                             &residual_index)
  {
    const Types::Index local_index =
      get_local_solution_index(residual_index, dependency_type);

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
    auto &feeval_variant = feeval_vector[(residual_index * max_dependency_types) +
                                         static_cast<Types::Index>(dependency_type)];

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
                          "Unexpected type in feeval_vector variant");
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
  FEEvaluationType                            *feeval_ptr,
  DiagonalType                                *diagonal_ptr,
  unsigned int                                 cell,
  Types::Index                                 global_variable_index,
  const std::function<void(VariableContainer &,
                           const dealii::Point<dim, SizeType> &,
                           const SizeType &)> &func,
  VectorType                                  &dst,
  const std::vector<VectorType *>             &src_subset)
{
  using DiagonalValueType = typename DiagonalType::value_type;

  // Grab the element volume
  const SizeType element_volume = element_volume_handler->get_element_volume(cell);

  // Helper function to submit the identity matrix
  auto submit_identity = [&](auto &_feeval_ptr, unsigned int dof_index)
  {
    for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
      {
        if constexpr (std::is_same_v<DiagonalValueType, SizeType> || dim == 1)
          {
            _feeval_ptr->submit_dof_value(SizeType(), j);
          }
        else
          {
            _feeval_ptr->submit_dof_value(DiagonalValueType(), j);
          }
      }

    // Set the i-th value to 1.0
    if constexpr (std::is_same_v<DiagonalValueType, SizeType> || dim == 1)
      {
        _feeval_ptr->submit_dof_value(dealii::make_vectorized_array<number>(1.0),
                                      dof_index);
      }
    else
      {
        DiagonalValueType one;
        for (unsigned int dimension = 0; dimension < dim; ++dimension)
          {
            one[dimension] = dealii::make_vectorized_array<number>(1.0);
          }
        _feeval_ptr->submit_dof_value(one, dof_index);
      }
  };
  // Reinit the cell for all the dependencies
  reinit(cell, global_variable_index);

  for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
    {
      // Submit an identity matrix for the change term
      submit_identity(feeval_ptr, i);

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

#include "core/variable_container.inst"

PRISMS_PF_END_NAMESPACE
