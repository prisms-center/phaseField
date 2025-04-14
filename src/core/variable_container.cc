// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/aligned_vector.h>
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

template <int dim, int degree, typename number>
variableContainer<dim, degree, number>::variableContainer(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
  const std::map<unsigned int, variableAttributes> &_subset_attributes,
  const std::map<std::pair<unsigned int, dependencyType>, unsigned int>
                  &_global_to_local_solution,
  const solveType &_solve_type)
  : subset_attributes(&_subset_attributes)
  , global_to_local_solution(&_global_to_local_solution)
  , solve_type(_solve_type)
{
  auto construct_map =
    [&](const std::map<unsigned int, std::map<dependencyType, fieldType>> &dependency_set)
  {
    for (const auto &[dependency_index, map] : dependency_set)
      {
        for (const auto &[dependency_type, field_type] : map)
          {
            if (field_type == fieldType::SCALAR)
              {
                feeval_map[dependency_index][dependency_type] =
                  std::make_unique<scalar_FEEval>(data, dependency_index);
              }
            else
              {
                feeval_map[dependency_index][dependency_type] =
                  std::make_unique<vector_FEEval>(data, dependency_index);
              }
          }
      }
  };

  // For explicit solves we have already flattened the dependencies
  if (solve_type == solveType::EXPLICIT_RHS || solve_type == solveType::POSTPROCESS)
    {
      construct_map(subset_attributes->begin()->second.dependency_set_RHS);
      return;
    }

  // TODO (landinjm): Add stuff for cononlinear solves

  // Loop through the variable attributes for nonexplicit solves
  Assert(subset_attributes->size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      construct_map(subset_attributes->begin()->second.dependency_set_LHS);
    }
  else
    {
      construct_map(subset_attributes->begin()->second.dependency_set_RHS);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::eval_local_operator(
  const std::function<void(variableContainer &, const dealii::Point<dim, size_type> &)>
                                              &func,
  std::vector<VectorType *>                   &dst,
  const std::vector<VectorType *>             &src,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      reinit_and_eval(src, cell);

      for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
        {
          // Set the quadrature point
          q_point = quad;

          // Grab the quadrature point location
          const dealii::Point<dim, size_type> q_point_loc = get_q_point_location();

          // Calculate the residuals
          func(*this, q_point_loc);
        }

      // Integrate and add to global vector dst
      integrate_and_distribute(dst);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::eval_local_operator(
  const std::function<void(variableContainer &, const dealii::Point<dim, size_type> &)>
                                              &func,
  VectorType                                  &dst,
  const std::vector<VectorType *>             &src,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      reinit_and_eval(src, cell);

      for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
        {
          // Set the quadrature point
          q_point = quad;

          // Grab the quadrature point location
          const dealii::Point<dim, size_type> q_point_loc = get_q_point_location();

          // Calculate the residuals
          func(*this, q_point_loc);
        }

      // Integrate and add to global vector dst
      integrate_and_distribute(dst);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::eval_local_operator(
  const std::function<void(variableContainer &, const dealii::Point<dim, size_type> &)>
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

      for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
        {
          // Set the quadrature point
          q_point = quad;

          // Grab the quadrature point location
          const dealii::Point<dim, size_type> q_point_loc = get_q_point_location();

          // Calculate the residuals
          func(*this, q_point_loc);
        }

      // Integrate and add to global vector dst
      integrate_and_distribute(dst);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::eval_local_diagonal(
  const std::function<void(variableContainer &, const dealii::Point<dim, size_type> &)>
                                              &func,
  VectorType                                  &dst,
  const std::vector<VectorType *>             &src_subset,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  Assert(subset_attributes->size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  const auto &global_var_index = subset_attributes->begin()->first;
  const auto &field_type       = subset_attributes->begin()->second.field_type;
  FEEval_exists(global_var_index, dependencyType::CHANGE);
  auto &feeval_variant = feeval_map.at(global_var_index).at(dependencyType::CHANGE);

  auto process_feeval = [&](auto &feeval_ptr, auto &diag_ptr)
  {
    using feeval_type      = std::decay_t<decltype(*feeval_ptr)>;
    using diag_type        = std::decay_t<decltype(*diag_ptr)>;
    using local_value_type = typename diag_type::value_type;

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        // Reinit the cell for all the dependencies
        reinit(cell, global_var_index);

        for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
          {
            // Submit an identity matrix for the change term
            for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
              {
                if constexpr (std::is_same_v<local_value_type, size_type> || dim == 1)
                  {
                    feeval_ptr->submit_dof_value(size_type(), j);
                  }
                else
                  {
                    feeval_ptr->submit_dof_value(local_value_type(), j);
                  }
              }

            // Set the i-th value to 1.0
            if constexpr (std::is_same_v<local_value_type, size_type> || dim == 1)
              {
                feeval_ptr->submit_dof_value(dealii::make_vectorized_array<number>(1.0),
                                             i);
              }
            else
              {
                local_value_type one;
                for (unsigned int dimension = 0; dimension < dim; ++dimension)
                  {
                    one[dimension] = dealii::make_vectorized_array<number>(1.0);
                  }
                feeval_ptr->submit_dof_value(one, i);
              }

            // Read plain dof values for non change src
            read_dof_values(src_subset);

            // Evaluate the dependencies based on the flags
            eval(global_var_index);

            for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
              {
                // Set the quadrature point
                q_point = quad;

                // Grab the quadrature point location
                const dealii::Point<dim, size_type> q_point_loc = get_q_point_location();

                // Calculate the residuals
                func(*this, q_point_loc);
              }

            // Integrate the diagonal
            integrate(global_var_index);
            (*diag_ptr)[i] = feeval_ptr->get_dof_value(i);
          }

        // Submit calculated diagonal values and distribute
        for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
          {
            if constexpr (std::is_same_v<local_value_type, size_type>)
              {
                feeval_ptr->submit_dof_value((*diag_ptr)[i], i);
              }
            else if constexpr (dim == 1)
              {
                feeval_ptr->submit_dof_value((*diag_ptr)[i][0], i);
              }
            else
              {
                feeval_ptr->submit_dof_value((*diag_ptr)[i], i);
              }
          }
        feeval_ptr->distribute_local_to_global(dst);
      }
  };

  if (field_type == fieldType::SCALAR)
    {
      scalar_FEEval *scalar_feeval_ptr = nullptr;

      if constexpr (dim == 1)
        {
          scalar_feeval_ptr = feeval_variant.get();
        }
      else
        {
          std::visit(
            [&scalar_feeval_ptr](auto &ptr)
            {
              using T = std::decay_t<decltype(ptr)>;
              if constexpr (std::is_same_v<T, std::unique_ptr<scalar_FEEval>>)
                {
                  scalar_feeval_ptr = ptr.get();
                }
              else
                {
                  Assert(false, dealii::ExcNotInitialized());
                }
            },
            feeval_variant);
        }

      n_dofs_per_cell              = scalar_feeval_ptr->dofs_per_cell;
      diagonal                     = std::make_unique<scalar_diag>(n_dofs_per_cell);
      scalar_diag *scalar_diag_ptr = nullptr;

      std::visit(
        [&scalar_diag_ptr](auto &ptr)
        {
          using T = std::decay_t<decltype(ptr)>;
          if constexpr (std::is_same_v<T, std::unique_ptr<scalar_diag>>)
            {
              scalar_diag_ptr = ptr.get();
            }
          else
            {
              Assert(false, dealii::ExcNotInitialized());
            }
        },
        diagonal);

      process_feeval(scalar_feeval_ptr, scalar_diag_ptr);
    }
  else if (field_type == fieldType::VECTOR)
    {
      vector_FEEval *vector_feeval_ptr = nullptr;

      if constexpr (dim == 1)
        {
          vector_feeval_ptr = feeval_variant.get();
        }
      else
        {
          std::visit(
            [&vector_feeval_ptr](auto &ptr)
            {
              using T = std::decay_t<decltype(ptr)>;
              if constexpr (std::is_same_v<T, std::unique_ptr<vector_FEEval>>)
                {
                  vector_feeval_ptr = ptr.get();
                }
              else
                {
                  Assert(false, dealii::ExcNotInitialized());
                }
            },
            feeval_variant);
        }

      n_dofs_per_cell              = vector_feeval_ptr->dofs_per_component;
      diagonal                     = std::make_unique<vector_diag>(n_dofs_per_cell);
      vector_diag *vector_diag_ptr = nullptr;

      std::visit(
        [&vector_diag_ptr](auto &ptr)
        {
          using T = std::decay_t<decltype(ptr)>;
          if constexpr (std::is_same_v<T, std::unique_ptr<vector_diag>>)
            {
              vector_diag_ptr = ptr.get();
            }
          else
            {
              Assert(false, dealii::ExcNotInitialized());
            }
        },
        diagonal);

      process_feeval(vector_feeval_ptr, vector_diag_ptr);
    }
  else
    {
      Assert(false, UnreachableCode());
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::FEEval_exists(
  [[maybe_unused]] const unsigned int   &dependency_index,
  [[maybe_unused]] const dependencyType &dependency_type) const
{
  Assert(feeval_map.contains(dependency_index),
         dealii::ExcMessage("The FEEvaluation object does not exist for global index = " +
                            std::to_string(dependency_index)));
  Assert(feeval_map.at(dependency_index).contains(dependency_type),
         dealii::ExcMessage("The FEEvaluation object with global index = " +
                            std::to_string(dependency_index) +
                            " does not exist for type = " + to_string(dependency_type)));
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::access_valid(
  [[maybe_unused]] const unsigned int                             &dependency_index,
  [[maybe_unused]] const dependencyType                           &dependency_type,
  [[maybe_unused]] const dealii::EvaluationFlags::EvaluationFlags &flag) const
{
  for ([[maybe_unused]] const auto &[index, variable] : *subset_attributes)
    {
      if (solve_type == solveType::NONEXPLICIT_LHS)
        {
          Assert(variable.eval_flag_set_LHS.contains(
                   std::make_pair(dependency_index, dependency_type)),
                 DependencyNotFound(dependency_index, to_string(dependency_type)));
        }
      else
        {
          Assert(variable.eval_flag_set_RHS.contains(
                   std::make_pair(dependency_index, dependency_type)),
                 DependencyNotFound(dependency_index, to_string(dependency_type)));
        }
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::submission_valid(
  [[maybe_unused]] const dependencyType &dependency_type) const
{
  Assert(dependency_type == NORMAL || solve_type == solveType::NONEXPLICIT_LHS,
         dealii::ExcMessage(
           "RHS residuals are only allowed to submit normal gradient terms."));
  Assert(dependency_type == CHANGE || solve_type != solveType::NONEXPLICIT_LHS,
         dealii::ExcMessage(
           "LHS residuals are only allowed to submit change gradient terms."));
}

template <int dim, int degree, typename number>
unsigned int
variableContainer<dim, degree, number>::get_n_q_points() const
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

template <int dim, int degree, typename number>
dealii::Point<dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_q_point_location() const
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
                [&](const auto &feeval_ptr) -> dealii::Point<dim, size_type>
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

  return dealii::Point<dim, size_type>();
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::reinit_and_eval(
  const std::vector<VectorType *> &src,
  unsigned int                     cell)
{
  // Reinit and eval values for the given dependency set. Note the dependency set includes
  // the variable we're evaluating, which may or may not be an actually dependency. For
  // this reason, I selectively read dofs and evaluate the flags.
  auto reinit_and_eval_map =
    [&](const std::map<std::pair<unsigned int, dependencyType>,
                       dealii::EvaluationFlags::EvaluationFlags>          &eval_flag_set,
        const std::map<unsigned int, std::map<dependencyType, fieldType>> &dependency_set)
  {
    for (const auto &[dependency_index, map] : dependency_set)
      {
        for (const auto &[dependency_type, field_type] : map)
          {
            if (dependency_type == dependencyType::CHANGE)
              {
                continue;
              }

            const auto &pair = std::make_pair(dependency_index, dependency_type);
            Assert(global_to_local_solution->contains(pair),
                   dealii::ExcMessage(
                     "The global to local mapping does not exists for global index = " +
                     std::to_string(dependency_index) +
                     "  and type = " + to_string(dependency_type)));
            FEEval_exists(dependency_index, dependency_type);
            auto &feeval_variant = feeval_map.at(dependency_index).at(dependency_type);

            auto process_feeval = [&](auto &feeval_ptr)
            {
              feeval_ptr->reinit(cell);
              if (!eval_flag_set.contains(pair))
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
              feeval_ptr->evaluate(eval_flag_set.at(pair));
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
                    static_assert(std::is_same_v<T, std::unique_ptr<scalar_FEEval>> ||
                                    std::is_same_v<T, std::unique_ptr<vector_FEEval>>,
                                  "Unexpected type in feeval_map variant");
                    process_feeval(ptr);
                  },
                  feeval_variant);
              }
          }
      }
  };

  if (solve_type == solveType::EXPLICIT_RHS || solve_type == solveType::POSTPROCESS)
    {
      const auto &attrs = subset_attributes->begin()->second;
      reinit_and_eval_map(attrs.eval_flag_set_RHS, attrs.dependency_set_RHS);
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
  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      reinit_and_eval_map(attrs.eval_flag_set_LHS, attrs.dependency_set_LHS);
    }
  else
    {
      reinit_and_eval_map(attrs.eval_flag_set_RHS, attrs.dependency_set_RHS);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::reinit_and_eval(const VectorType &src,
                                                        unsigned int      cell)
{
  // Reinit and eval values for the given dependency set. Note the dependency set includes
  // the variable we're evaluating, which may or may not be an actually dependency. For
  // this reason, I selectively read dofs and evaluate the flags.
  auto reinit_and_eval_map =
    [&](const std::map<std::pair<unsigned int, dependencyType>,
                       dealii::EvaluationFlags::EvaluationFlags>          &eval_flag_set,
        const std::map<unsigned int, std::map<dependencyType, fieldType>> &dependency_set)
  {
    for (const auto &[dependency_index, map] : dependency_set)
      {
        for (const auto &[dependency_type, field_type] : map)
          {
            const auto &pair = std::make_pair(dependency_index, dependency_type);
            Assert(global_to_local_solution->contains(pair),
                   dealii::ExcMessage(
                     "The global to local mapping does not exists for global index = " +
                     std::to_string(dependency_index) +
                     "  and type = " + to_string(dependency_type)));
            FEEval_exists(dependency_index, dependency_type);
            auto &feeval_variant = feeval_map.at(dependency_index).at(dependency_type);

            auto process_feeval = [&](auto &feeval_ptr)
            {
              feeval_ptr->reinit(cell);
              feeval_ptr->read_dof_values_plain(src);
              feeval_ptr->evaluate(eval_flag_set.at(pair));
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
                    static_assert(std::is_same_v<T, std::unique_ptr<scalar_FEEval>> ||
                                    std::is_same_v<T, std::unique_ptr<vector_FEEval>>,
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
  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      reinit_and_eval_map(attrs.eval_flag_set_LHS, attrs.dependency_set_LHS);
      return;
    }
  Assert(false,
         dealii::ExcMessage(
           "reinit_and_eval(src) should only be called for LHS evaluations"));
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::reinit(unsigned int        cell,
                                               const unsigned int &global_variable_index)
{
  auto reinit_map =
    [&](const std::map<std::pair<unsigned int, dependencyType>,
                       dealii::EvaluationFlags::EvaluationFlags> &eval_flag_set)
  {
    for (const auto &[pair, flags] : eval_flag_set)
      {
        const unsigned int   &dependency_index = pair.first;
        const dependencyType &dependency_type  = pair.second;

        Assert(subset_attributes->contains(dependency_index),
               dealii::ExcMessage(
                 "The subset attribute entry does not exists for global index = " +
                 std::to_string(dependency_index)));
        FEEval_exists(dependency_index, dependency_type);
        auto &feeval_variant = feeval_map.at(dependency_index).at(dependency_type);

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
                static_assert(std::is_same_v<T, std::unique_ptr<scalar_FEEval>> ||
                                std::is_same_v<T, std::unique_ptr<vector_FEEval>>,
                              "Unexpected type in feeval_map variant");
                process_feeval(ptr);
              },
              feeval_variant);
          }
      }
  };

  Assert(subset_attributes->contains(global_variable_index),
         dealii::ExcMessage(
           "The subset attribute entry does not exists for global index = " +
           std::to_string(global_variable_index)));
  const auto &attrs = subset_attributes->at(global_variable_index);
  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      reinit_map(attrs.eval_flag_set_LHS);
    }
  else
    {
      reinit_map(attrs.eval_flag_set_RHS);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::read_dof_values(
  const std::vector<VectorType *> &src)
{
  auto reinit_and_eval_map =
    [&](const std::map<std::pair<unsigned int, dependencyType>,
                       dealii::EvaluationFlags::EvaluationFlags>          &eval_flag_set,
        const std::map<unsigned int, std::map<dependencyType, fieldType>> &dependency_set)
  {
    for (const auto &[dependency_index, map] : dependency_set)
      {
        for (const auto &[dependency_type, field_type] : map)
          {
            if (dependency_type == dependencyType::CHANGE)
              {
                continue;
              }

            const auto &pair = std::make_pair(dependency_index, dependency_type);
            Assert(global_to_local_solution->contains(pair),
                   dealii::ExcMessage(
                     "The global to local mapping does not exists for global index = " +
                     std::to_string(dependency_index) +
                     "  and type = " + to_string(dependency_type)));
            FEEval_exists(dependency_index, dependency_type);
            auto &feeval_variant = feeval_map.at(dependency_index).at(dependency_type);

            auto process_feeval = [&](auto &feeval_ptr)
            {
              if (eval_flag_set.contains(pair))
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
                    static_assert(std::is_same_v<T, std::unique_ptr<scalar_FEEval>> ||
                                    std::is_same_v<T, std::unique_ptr<vector_FEEval>>,
                                  "Unexpected type in feeval_map variant");
                    process_feeval(ptr);
                  },
                  feeval_variant);
              }
          }
      }
  };

  if (solve_type != solveType::NONEXPLICIT_LHS)
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
  reinit_and_eval_map(attrs.eval_flag_set_LHS, attrs.dependency_set_LHS);
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::eval(const unsigned int &global_variable_index)
{
  auto eval_map =
    [&](const std::map<std::pair<unsigned int, dependencyType>,
                       dealii::EvaluationFlags::EvaluationFlags> &eval_flag_set)
  {
    for (const auto &[pair, flags] : eval_flag_set)
      {
        const unsigned int   &dependency_index = pair.first;
        const dependencyType &dependency_type  = pair.second;

        Assert(subset_attributes->contains(dependency_index),
               dealii::ExcMessage(
                 "The subset attribute entry does not exists for global index = " +
                 std::to_string(dependency_index)));
        FEEval_exists(dependency_index, dependency_type);
        auto &feeval_variant = feeval_map.at(dependency_index).at(dependency_type);

        auto process_feeval = [&](auto &feeval_ptr)
        {
          feeval_ptr->evaluate(flags);
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
                static_assert(std::is_same_v<T, std::unique_ptr<scalar_FEEval>> ||
                                std::is_same_v<T, std::unique_ptr<vector_FEEval>>,
                              "Unexpected type in feeval_map variant");
                process_feeval(ptr);
              },
              feeval_variant);
          }
      }
  };

  Assert(subset_attributes->contains(global_variable_index),
         dealii::ExcMessage(
           "The subset attribute entry does not exists for global index = " +
           std::to_string(global_variable_index)));
  const auto &attrs = subset_attributes->at(global_variable_index);
  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      eval_map(attrs.eval_flag_set_LHS);
    }
  else
    {
      eval_map(attrs.eval_flag_set_RHS);
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::integrate(
  const unsigned int &global_variable_index)
{
  Assert(subset_attributes->contains(global_variable_index),
         dealii::ExcMessage(
           "The subset attribute entry does not exists for global index = " +
           std::to_string(global_variable_index)));
  const auto &variable = subset_attributes->at(global_variable_index);

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      FEEval_exists(global_variable_index, dependencyType::CHANGE);
      auto &feeval_variant =
        feeval_map.at(global_variable_index).at(dependencyType::CHANGE);

      auto process_feeval = [&](auto &feeval_ptr)
      {
        feeval_ptr->integrate(variable.eval_flags_residual_LHS);
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
              static_assert(std::is_same_v<T, std::unique_ptr<scalar_FEEval>> ||
                              std::is_same_v<T, std::unique_ptr<vector_FEEval>>,
                            "Unexpected type in feeval_map variant");
              process_feeval(ptr);
            },
            feeval_variant);
        }
      return;
    }

  AssertThrow(false,
              dealii::ExcMessage(
                "Integrate called for a solve type that is not NONEXPLICIT_LHS."));
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::integrate_and_distribute(
  std::vector<VectorType *> &dst)
{
  auto integrate_and_distribute_map =
    [&](const dealii::EvaluationFlags::EvaluationFlags &residual_flag_set,
        const dependencyType                           &dependency_type,
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
    FEEval_exists(residual_index, dependency_type);
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
            static_assert(std::is_same_v<T, std::unique_ptr<scalar_FEEval>> ||
                            std::is_same_v<T, std::unique_ptr<vector_FEEval>>,
                          "Unexpected type in feeval_map variant");
            process_feeval(ptr);
          },
          feeval_variant);
      }
  };

  for (const auto &[index, variable] : *subset_attributes)
    {
      if (solve_type == solveType::NONEXPLICIT_LHS)
        {
          integrate_and_distribute_map(variable.eval_flags_residual_LHS,
                                       dependencyType::CHANGE,
                                       index);
        }
      else
        {
          integrate_and_distribute_map(variable.eval_flags_residual_RHS,
                                       dependencyType::NORMAL,
                                       index);
        }
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::integrate_and_distribute(VectorType &dst)
{
  auto integrate_and_distribute_map =
    [&](const dealii::EvaluationFlags::EvaluationFlags &residual_flag_set,
        const dependencyType                           &dependency_type,
        const unsigned int                             &residual_index)
  {
    Assert(subset_attributes->contains(residual_index),
           dealii::ExcMessage(
             "The subset attribute entry does not exists for global index = " +
             std::to_string(residual_index)));
    FEEval_exists(residual_index, dependency_type);
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
            static_assert(std::is_same_v<T, std::unique_ptr<scalar_FEEval>> ||
                            std::is_same_v<T, std::unique_ptr<vector_FEEval>>,
                          "Unexpected type in feeval_map variant");
            process_feeval(ptr);
          },
          feeval_variant);
      }
  };

  Assert(subset_attributes->size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  if (solve_type == solveType::NONEXPLICIT_LHS)
    {
      integrate_and_distribute_map(
        subset_attributes->begin()->second.eval_flags_residual_LHS,
        dependencyType::CHANGE,
        subset_attributes->begin()->first);
    }
  else
    {
      integrate_and_distribute_map(
        subset_attributes->begin()->second.eval_flags_residual_RHS,
        dependencyType::NORMAL,
        subset_attributes->begin()->first);
    }
}

template <int dim, int degree, typename number>
typename variableContainer<dim, degree, number>::size_type
variableContainer<dim, degree, number>::get_scalar_value(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::values);
  FEEval_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      return feeval_map.at(global_variable_index).at(dependency_type)->get_value(q_point);
    }
  else
    {
      return std::visit(
        [&](const auto &feeval_ptr) -> size_type
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, scalar_FEEval>)
            {
              return feeval_ptr->get_value(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected scalar_FEEval but got vector_FEEval."));
              return size_type();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_scalar_gradient(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  FEEval_exists(global_variable_index, dependency_type);
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
        [&](const auto &feeval_ptr) -> dealii::Tensor<1, dim, size_type>
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, scalar_FEEval>)
            {
              return feeval_ptr->get_gradient(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected scalar_FEEval but got vector_FEEval."));
              return dealii::Tensor<1, dim, size_type>();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_scalar_hessian(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  FEEval_exists(global_variable_index, dependency_type);
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
        [&](const auto &feeval_ptr) -> dealii::Tensor<2, dim, size_type>
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, scalar_FEEval>)
            {
              return feeval_ptr->get_hessian(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected scalar_FEEval but got vector_FEEval."));
              return dealii::Tensor<2, dim, size_type>();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_scalar_hessian_diagonal(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  FEEval_exists(global_variable_index, dependency_type);
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
        [&](const auto &feeval_ptr) -> dealii::Tensor<1, dim, size_type>
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, scalar_FEEval>)
            {
              return feeval_ptr->get_hessian_diagonal(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected scalar_FEEval but got vector_FEEval."));
              return dealii::Tensor<1, dim, size_type>();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
typename variableContainer<dim, degree, number>::size_type
variableContainer<dim, degree, number>::get_scalar_laplacian(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  FEEval_exists(global_variable_index, dependency_type);
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
        [&](const auto &feeval_ptr) -> size_type
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, scalar_FEEval>)
            {
              return feeval_ptr->get_laplacian(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected scalar_FEEval but got vector_FEEval."));
              return size_type();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_value(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::values);
  FEEval_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      dealii::Tensor<1, dim, size_type> wrapper;
      wrapper[0] =
        feeval_map.at(global_variable_index).at(dependency_type)->get_value(q_point);
      return wrapper;
    }
  else
    {
      return std::visit(
        [&](const auto &feeval_ptr) -> dealii::Tensor<1, dim, size_type>
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, vector_FEEval>)
            {
              return feeval_ptr->get_value(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected vector_FEEval but got scalar_FEEval."));
              return dealii::Tensor<1, dim, size_type>();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_gradient(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  FEEval_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      dealii::Tensor<2, dim, size_type> wrapper;
      wrapper[0] =
        feeval_map.at(global_variable_index).at(dependency_type)->get_gradient(q_point);
      return wrapper;
    }
  else
    {
      return std::visit(
        [&](const auto &feeval_ptr) -> dealii::Tensor<2, dim, size_type>
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, vector_FEEval>)
            {
              return feeval_ptr->get_gradient(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected vector_FEEval but got scalar_FEEval."));
              return dealii::Tensor<2, dim, size_type>();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<3, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_hessian(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  FEEval_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      dealii::Tensor<3, dim, size_type> wrapper;
      wrapper[0] =
        feeval_map.at(global_variable_index).at(dependency_type)->get_hessian(q_point);
      return wrapper;
    }
  else
    {
      return std::visit(
        [&](const auto &feeval_ptr) -> dealii::Tensor<3, dim, size_type>
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, vector_FEEval>)
            {
              return feeval_ptr->get_hessian(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected vector_FEEval but got scalar_FEEval."));
              return dealii::Tensor<3, dim, size_type>();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_hessian_diagonal(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  FEEval_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      dealii::Tensor<2, dim, size_type> wrapper;
      wrapper[0] = feeval_map.at(global_variable_index)
                     .at(dependency_type)
                     ->get_hessian_diagonal(q_point);
      return wrapper;
    }
  else
    {
      return std::visit(
        [&](const auto &feeval_ptr) -> dealii::Tensor<2, dim, size_type>
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, vector_FEEval>)
            {
              return feeval_ptr->get_hessian_diagonal(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected vector_FEEval but got scalar_FEEval."));
              return dealii::Tensor<2, dim, size_type>();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<1, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_laplacian(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::hessians);
  FEEval_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      dealii::Tensor<1, dim, size_type> wrapper;
      wrapper[0] =
        feeval_map.at(global_variable_index).at(dependency_type)->get_laplacian(q_point);
      return wrapper;
    }
  else
    {
      return std::visit(
        [&](const auto &feeval_ptr) -> dealii::Tensor<1, dim, size_type>
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, vector_FEEval>)
            {
              return feeval_ptr->get_laplacian(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected vector_FEEval but got scalar_FEEval."));
              return dealii::Tensor<1, dim, size_type>();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
typename variableContainer<dim, degree, number>::size_type
variableContainer<dim, degree, number>::get_vector_divergence(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  FEEval_exists(global_variable_index, dependency_type);
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
        [&](const auto &feeval_ptr) -> size_type
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, vector_FEEval>)
            {
              return feeval_ptr->get_divergence(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected vector_FEEval but got scalar_FEEval."));
              return size_type();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<2, dim, typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_symmetric_gradient(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  FEEval_exists(global_variable_index, dependency_type);
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
        [&](const auto &feeval_ptr) -> dealii::Tensor<2, dim, size_type>
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, vector_FEEval>)
            {
              return feeval_ptr->get_symmetric_gradient(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected vector_FEEval but got scalar_FEEval."));
              return dealii::Tensor<2, dim, size_type>();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
dealii::Tensor<1,
               (dim == 2 ? 1 : dim),
               typename variableContainer<dim, degree, number>::size_type>
variableContainer<dim, degree, number>::get_vector_curl(
  unsigned int   global_variable_index,
  dependencyType dependency_type) const
{
#ifdef DEBUG
  access_valid(global_variable_index,
               dependency_type,
               dealii::EvaluationFlags::EvaluationFlags::gradients);
  FEEval_exists(global_variable_index, dependency_type);
#endif

  if constexpr (dim == 1)
    {
      Assert(false, dealii::ExcMessage("Curl is nonsensical for 1D."));
      return dealii::Tensor<1, (dim == 2 ? 1 : dim), size_type> {};
    }
  else
    {
      // Return the value directly for dim > 1
      return std::visit(
        [&](const auto &feeval_ptr) -> dealii::Tensor<1, (dim == 2 ? 1 : dim), size_type>
        {
          using FEEvalType = std::decay_t<decltype(*feeval_ptr)>;
          if constexpr (std::is_same_v<FEEvalType, vector_FEEval>)
            {
              return feeval_ptr->get_curl(q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected vector_FEEval but got scalar_FEEval."));
              return dealii::Tensor<1, (dim == 2 ? 1 : dim), size_type>();
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_scalar_value_term(
  const unsigned int   &global_variable_index,
  const size_type      &val,
  const dependencyType &dependency_type)
{
#ifdef DEBUG
  submission_valid(dependency_type);
  FEEval_exists(global_variable_index, dependency_type);
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
          if constexpr (std::is_same_v<FEEvalType, scalar_FEEval>)
            {
              feeval_ptr->submit_value(val, q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected scalar_FEEval but got vector_FEEval."));
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_scalar_gradient_term(
  const unsigned int                      &global_variable_index,
  const dealii::Tensor<1, dim, size_type> &grad,
  const dependencyType                    &dependency_type)
{
#ifdef DEBUG
  submission_valid(dependency_type);
  FEEval_exists(global_variable_index, dependency_type);
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
          if constexpr (std::is_same_v<FEEvalType, scalar_FEEval>)
            {
              feeval_ptr->submit_gradient(grad, q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected scalar_FEEval but got vector_FEEval."));
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_vector_value_term(
  const unsigned int                      &global_variable_index,
  const dealii::Tensor<1, dim, size_type> &val,
  const dependencyType                    &dependency_type)
{
#ifdef DEBUG
  submission_valid(dependency_type);
  FEEval_exists(global_variable_index, dependency_type);
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
          if constexpr (std::is_same_v<FEEvalType, vector_FEEval>)
            {
              feeval_ptr->submit_value(val, q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected vector_FEEval but got scalar_FEEval."));
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

template <int dim, int degree, typename number>
void
variableContainer<dim, degree, number>::set_vector_gradient_term(
  const unsigned int                      &global_variable_index,
  const dealii::Tensor<2, dim, size_type> &grad,
  const dependencyType                    &dependency_type)
{
#ifdef DEBUG
  submission_valid(dependency_type);
  FEEval_exists(global_variable_index, dependency_type);
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
          if constexpr (std::is_same_v<FEEvalType, vector_FEEval>)
            {
              feeval_ptr->submit_gradient(grad, q_point);
            }
          else
            {
              Assert(false,
                     dealii::ExcMessage("Expected vector_FEEval but got scalar_FEEval."));
            }
        },
        feeval_map.at(global_variable_index).at(dependency_type));
    }
}

INSTANTIATE_TRI_TEMPLATE(variableContainer)

PRISMS_PF_END_NAMESPACE
