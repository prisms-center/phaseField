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
  const std::map<std::pair<unsigned int, DependencyType>, unsigned int>
                  &_global_to_local_solution,
  const SolveType &_solve_type,
  bool             use_local_mapping)
  : subset_attributes(&_subset_attributes)
  , global_to_local_solution(&_global_to_local_solution)
  , solve_type(_solve_type)
{
  auto construct_map =
    [&](const std::map<unsigned int, std::map<DependencyType, FieldType>> &dependency_set)
  {
    for (const auto &[dependency_index, map] : dependency_set)
      {
        for (const auto &[dependency_type, field_type] : map)
          {
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
                    feeval_map[dependency_index][dependency_type] =
                      std::make_unique<ScalarFEEvaluation>(
                        data,
                        global_to_local_solution->at(
                          std::make_pair(dependency_index, dependency_type)));
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
                    feeval_map[dependency_index][dependency_type] =
                      std::make_unique<VectorFEEvaluation>(
                        data,
                        global_to_local_solution->at(
                          std::make_pair(dependency_index, dependency_type)));
                  }
              }
          }
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
  auto &feeval_variant = feeval_map.at(global_var_index).at(DependencyType::Change);

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
  Assert(feeval_map.contains(dependency_index),
         dealii::ExcMessage("The FEEvaluation object does not exist for global index = " +
                            std::to_string(dependency_index)));
  Assert(feeval_map.at(dependency_index).contains(dependency_type),
         dealii::ExcMessage("The FEEvaluation object with global index = " +
                            std::to_string(dependency_index) +
                            " does not exist for type = " + to_string(dependency_type)));
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
  for (const auto &[dependency_index, dependency_map] : feeval_map)
    {
      for (const auto &[dependency_type, feeval_variant] : dependency_map)
        {
          if constexpr (dim == 1)
            {
              return feeval_variant->n_q_points;
            }
          else
            {
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
  for (const auto &[dependency_index, dependency_map] : feeval_map)
    {
      for (const auto &[dependency_type, feeval_variant] : dependency_map)
        {
          if constexpr (dim == 1)
            {
              return feeval_variant->quadrature_point(q_point);
            }
          else
            {
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
        const std::map<unsigned int, std::map<DependencyType, FieldType>> &dependency_set)
  {
    for (const auto &[dependency_index, map] : dependency_set)
      {
        for (const auto &[dependency_type, field_type] : map)
          {
            if (dependency_type == DependencyType::Change)
              {
                continue;
              }

            const auto &pair = std::make_pair(dependency_index, dependency_type);
            Assert(global_to_local_solution->contains(pair),
                   dealii::ExcMessage(
                     "The global to local mapping does not exists for global index = " +
                     std::to_string(dependency_index) +
                     "  and type = " + to_string(dependency_type)));
            feevaluation_exists(dependency_index, dependency_type);
            auto &feeval_variant = feeval_map.at(dependency_index).at(dependency_type);

            auto process_feeval = [&](auto &feeval_ptr)
            {
              feeval_ptr->reinit(cell);
              if (eval_flag_set[dependency_index]
                               [static_cast<Types::Index>(dependency_type)] ==
                  dealii::EvaluationFlags::EvaluationFlags::nothing)
                {
                  return;
                }
              const unsigned int &local_index = global_to_local_solution->at(pair);
              Assert(src.size() > local_index,
                     dealii::ExcMessage(
                       "The provided src vector's size is below the given local "
                       "index = " +
                       std::to_string(local_index) +
                       " for global index = " + std::to_string(dependency_index) +
                       "  and type = " + to_string(dependency_type)));
              feeval_ptr->read_dof_values_plain(*(src.at(local_index)));
              feeval_ptr->evaluate(
                eval_flag_set[dependency_index]
                             [static_cast<Types::Index>(dependency_type)]);
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
          }
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
        const std::map<unsigned int, std::map<DependencyType, FieldType>> &dependency_set)
  {
    for (const auto &[dependency_index, map] : dependency_set)
      {
        for (const auto &[dependency_type, field_type] : map)
          {
            // TODO (landinjm): This can be drastically simplified because all we're doing
            // in reinit-ing and eval-ing the change solution.
            if (dependency_type != DependencyType::Change)
              {
                continue;
              }

            feevaluation_exists(dependency_index, dependency_type);
            auto &feeval_variant = feeval_map.at(dependency_index).at(dependency_type);

            auto process_feeval = [&](auto &feeval_ptr)
            {
              feeval_ptr->reinit(cell);
              feeval_ptr->read_dof_values_plain(src);
              feeval_ptr->evaluate(
                eval_flag_set[dependency_index]
                             [static_cast<Types::Index>(dependency_type)]);
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
          }
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
            auto &feeval_variant = feeval_map.at(dependency_index)
                                     .at(static_cast<DependencyType>(dependency_type));

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
        const std::map<unsigned int, std::map<DependencyType, FieldType>> &dependency_set)
  {
    for (const auto &[dependency_index, map] : dependency_set)
      {
        for (const auto &[dependency_type, field_type] : map)
          {
            if (dependency_type == DependencyType::Change)
              {
                continue;
              }

            const auto &pair = std::make_pair(dependency_index, dependency_type);
            Assert(global_to_local_solution->contains(pair),
                   dealii::ExcMessage(
                     "The global to local mapping does not exists for global index = " +
                     std::to_string(dependency_index) +
                     "  and type = " + to_string(dependency_type)));
            feevaluation_exists(dependency_index, dependency_type);
            auto &feeval_variant = feeval_map.at(dependency_index).at(dependency_type);

            auto process_feeval = [&](auto &feeval_ptr)
            {
              if (eval_flag_set[dependency_index]
                               [static_cast<Types::Index>(dependency_type)] !=
                  dealii::EvaluationFlags::EvaluationFlags::nothing)
                {
                  const unsigned int &local_index = global_to_local_solution->at(pair);
                  Assert(src.size() > local_index,
                         dealii::ExcMessage(
                           "The provided src vector's size is below the given local "
                           "index = " +
                           std::to_string(local_index) +
                           " for global index = " + std::to_string(dependency_index) +
                           "  and type = " + to_string(dependency_type)));
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
          }
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
            auto &feeval_variant =
              feeval_map.at(index).at(static_cast<DependencyType>(dep_index));

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
        feeval_map.at(global_variable_index).at(DependencyType::Change);

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
    Assert(
      global_to_local_solution->contains(std::make_pair(residual_index, dependency_type)),
      dealii::ExcMessage(
        "The global to local mapping does not exists for global index = " +
        std::to_string(residual_index) + "  and type = " + to_string(dependency_type)));
    const unsigned int &local_index =
      global_to_local_solution->at(std::make_pair(residual_index, dependency_type));
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
    auto &feeval_variant = feeval_map.at(residual_index).at(dependency_type);

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
    auto &feeval_variant = feeval_map.at(residual_index).at(dependency_type);

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
VariableContainer<dim, degree, number>::get_scalar_value(
  unsigned int   global_variable_index,
  DependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::values);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      return feeval_map.at(global_variable_index).at(dependency_type)->get_value(q_point);
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
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Tensor<1, dim, typename VariableContainer<dim, degree, number>::SizeType>
VariableContainer<dim, degree, number>::get_scalar_gradient(
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
      return feeval_map.at(global_variable_index)
        .at(dependency_type)
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
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Tensor<2, dim, typename VariableContainer<dim, degree, number>::SizeType>
VariableContainer<dim, degree, number>::get_scalar_hessian(
  unsigned int   global_variable_index,
  DependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      return feeval_map.at(global_variable_index)
        .at(dependency_type)
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
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Tensor<1, dim, typename VariableContainer<dim, degree, number>::SizeType>
VariableContainer<dim, degree, number>::get_scalar_hessian_diagonal(
  unsigned int   global_variable_index,
  DependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      return feeval_map.at(global_variable_index)
        .at(dependency_type)
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
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
typename VariableContainer<dim, degree, number>::SizeType
VariableContainer<dim, degree, number>::get_scalar_laplacian(
  unsigned int   global_variable_index,
  DependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      return feeval_map.at(global_variable_index)
        .at(dependency_type)
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
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Tensor<1, dim, typename VariableContainer<dim, degree, number>::SizeType>
VariableContainer<dim, degree, number>::get_vector_value(
  unsigned int   global_variable_index,
  DependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::values);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      dealii::Tensor<1, dim, SizeType> wrapper;
      wrapper[0] =
        feeval_map.at(global_variable_index).at(dependency_type)->get_value(q_point);
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
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Tensor<2, dim, typename VariableContainer<dim, degree, number>::SizeType>
VariableContainer<dim, degree, number>::get_vector_gradient(
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
      dealii::Tensor<2, dim, SizeType> wrapper;
      wrapper[0] =
        feeval_map.at(global_variable_index).at(dependency_type)->get_gradient(q_point);
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
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Tensor<3, dim, typename VariableContainer<dim, degree, number>::SizeType>
VariableContainer<dim, degree, number>::get_vector_hessian(
  unsigned int   global_variable_index,
  DependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      dealii::Tensor<3, dim, SizeType> wrapper;
      wrapper[0] =
        feeval_map.at(global_variable_index).at(dependency_type)->get_hessian(q_point);
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
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Tensor<2, dim, typename VariableContainer<dim, degree, number>::SizeType>
VariableContainer<dim, degree, number>::get_vector_hessian_diagonal(
  unsigned int   global_variable_index,
  DependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      dealii::Tensor<2, dim, SizeType> wrapper;
      wrapper[0] = feeval_map.at(global_variable_index)
                     .at(dependency_type)
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
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Tensor<1, dim, typename VariableContainer<dim, degree, number>::SizeType>
VariableContainer<dim, degree, number>::get_vector_laplacian(
  unsigned int   global_variable_index,
  DependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      dealii::Tensor<1, dim, SizeType> wrapper;
      wrapper[0] =
        feeval_map.at(global_variable_index).at(dependency_type)->get_laplacian(q_point);
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
        feeval_map.at(global_variable_index).at(dependency_type));
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
      return feeval_map.at(global_variable_index)
        .at(dependency_type)
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
        feeval_map.at(global_variable_index).at(dependency_type));
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
      return feeval_map.at(global_variable_index)
        .at(dependency_type)
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
        feeval_map.at(global_variable_index).at(dependency_type));
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
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::set_scalar_value_term(
  const unsigned int   &global_variable_index,
  const SizeType       &val,
  const DependencyType &dependency_type)
{
#ifdef DEBUG
  submission_valid(dependency_type);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      return feeval_map.at(global_variable_index)
        .at(dependency_type)
        ->submit_value(val, q_point);
    }
  else
    {
      std::visit(
        [&](auto &feeval_ptr)
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
            {
              feeval_ptr->submit_value(val, q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage(
                       "Expected ScalarFEEvaluation but got VectorFEEvaluation."));
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::set_scalar_gradient_term(
  const unsigned int                     &global_variable_index,
  const dealii::Tensor<1, dim, SizeType> &grad,
  const DependencyType                   &dependency_type)
{
#ifdef DEBUG
  submission_valid(dependency_type);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      return feeval_map.at(global_variable_index)
        .at(dependency_type)
        ->submit_gradient(grad, q_point);
    }
  else
    {
      std::visit(
        [&](auto &feeval_ptr)
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, ScalarFEEvaluation>)
            {
              feeval_ptr->submit_gradient(grad, q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage(
                       "Expected ScalarFEEvaluation but got VectorFEEvaluation."));
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::set_vector_value_term(
  const unsigned int                     &global_variable_index,
  const dealii::Tensor<1, dim, SizeType> &val,
  const DependencyType                   &dependency_type)
{
#ifdef DEBUG
  submission_valid(dependency_type);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      return feeval_map.at(global_variable_index)
        .at(dependency_type)
        ->submit_value(val, q_point);
    }
  else
    {
      std::visit(
        [&](auto &feeval_ptr)
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
            {
              feeval_ptr->submit_value(val, q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage(
                       "Expected VectorFEEvaluation but got ScalarFEEvaluation."));
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
VariableContainer<dim, degree, number>::set_vector_gradient_term(
  const unsigned int                     &global_variable_index,
  const dealii::Tensor<2, dim, SizeType> &grad,
  const DependencyType                   &dependency_type)
{
#ifdef DEBUG
  submission_valid(dependency_type);
  feevaluation_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      return feeval_map.at(global_variable_index)
        .at(dependency_type)
        ->submit_gradient(grad, q_point);
    }
  else
    {
      std::visit(
        [&](auto &feeval_ptr)
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, VectorFEEvaluation>)
            {
              feeval_ptr->submit_gradient(grad, q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage(
                       "Expected VectorFEEvaluation but got ScalarFEEvaluation."));
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

INSTANTIATE_TRI_TEMPLATE(VariableContainer)

PRISMS_PF_END_NAMESPACE
