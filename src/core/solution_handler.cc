// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <array>
#include <iterator>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
SolutionHandler<dim, number>::SolutionHandler(
  const std::map<unsigned int, VariableAttributes> &_attributes_list,
  const MGInfo<dim>                                &_mg_info)
  : attributes_list(&_attributes_list)
  , mg_info(&_mg_info)
{
  // If we don't have multigrid, we can return early
  if (!_mg_info.has_multigrid())
    {
      return;
    }
  has_multigrid = true;

  global_min_level = _mg_info.get_mg_min_level();
  mg_solution_set.resize(_mg_info.get_mg_depth());
  for (unsigned int level = 0; level < _mg_info.get_mg_depth(); level++)
    {
      mg_solution_set[level].resize(_mg_info.get_mg_breadth(level));
    }
}

template <unsigned int dim, typename number>
std::map<unsigned int, typename SolutionHandler<dim, number>::VectorType *>
SolutionHandler<dim, number>::get_solution_vector() const
{
  std::map<unsigned int, VectorType *> temp;
  for (const auto &[pair, vector] : solution_set)
    {
      if (pair.second != DependencyType::Normal)
        {
          continue;
        }

      Assert(vector != nullptr, dealii::ExcNotInitialized());

      temp.emplace(pair.first, vector.get());
    }

  return temp;
}

template <unsigned int dim, typename number>
typename SolutionHandler<dim, number>::VectorType *
SolutionHandler<dim, number>::get_solution_vector(unsigned int   index,
                                                  DependencyType dependency_type) const
{
  const auto pair = std::make_pair(index, dependency_type);

  Assert(solution_set.contains(pair),
         dealii::ExcMessage(
           "There is no solution vector for the given index = " + std::to_string(index) +
           " and type = " + to_string(dependency_type)));
  Assert(solution_set.at(pair) != nullptr, dealii::ExcNotInitialized());

  return solution_set.at(pair).get();
}

template <unsigned int dim, typename number>
std::map<unsigned int, typename SolutionHandler<dim, number>::VectorType *>
SolutionHandler<dim, number>::get_new_solution_vector() const
{
  std::map<unsigned int, VectorType *> temp;
  for (const auto &[index, vector] : new_solution_set)
    {
      Assert(vector != nullptr, dealii::ExcNotInitialized());

      temp.emplace(index, vector.get());
    }

  return temp;
}

template <unsigned int dim, typename number>
typename SolutionHandler<dim, number>::VectorType *
SolutionHandler<dim, number>::get_new_solution_vector(unsigned int index) const
{
  Assert(new_solution_set.contains(index),
         dealii::ExcMessage("There is no new solution vector for the given index = " +
                            std::to_string(index)));
  Assert(new_solution_set.at(index) != nullptr, dealii::ExcNotInitialized());

  return new_solution_set.at(index).get();
}

template <unsigned int dim, typename number>
std::vector<typename SolutionHandler<dim, number>::MGVectorType *>
SolutionHandler<dim, number>::get_mg_solution_vector(unsigned int level) const
{
  Assert(has_multigrid, dealii::ExcNotInitialized());

  // Convert absolute level to a relative level
  const unsigned int relative_level = level - global_min_level;
  Assert(relative_level < mg_solution_set.size(),
         dealii::ExcMessage("The mg solution set does not contain level = " +
                            std::to_string(level)));
  std::vector<MGVectorType *> temp;
  temp.reserve(mg_solution_set[relative_level].size());

  std::transform(mg_solution_set[relative_level].begin(),
                 mg_solution_set[relative_level].end(),
                 std::back_inserter(temp),
                 [](const std::unique_ptr<MGVectorType> &vector)
                 {
                   return vector.get();
                 });
  return temp;
}

template <unsigned int dim, typename number>
typename SolutionHandler<dim, number>::MGVectorType *
SolutionHandler<dim, number>::get_mg_solution_vector(unsigned int level,
                                                     unsigned int index) const
{
  Assert(has_multigrid, dealii::ExcNotInitialized());

  // Convert absolute level to a relative level
  const unsigned int relative_level = level - global_min_level;
  Assert(relative_level < mg_solution_set.size(),
         dealii::ExcMessage("The mg solution set does not contain level = " +
                            std::to_string(level)));
  [[maybe_unused]] const unsigned int global_index =
    mg_info->get_global_index(index, relative_level);
  Assert(index < mg_solution_set[relative_level].size(),
         dealii::ExcMessage(
           "The mg solution at the given level does not contain index = " +
           std::to_string(global_index)));
  return mg_solution_set[relative_level][index].get();
}

template <unsigned int dim, typename number>
void
SolutionHandler<dim, number>::init(MatrixfreeHandler<dim, number> &matrix_free_handler)
{
  // Create all entries
  for (const auto &[index, variable] : *attributes_list)
    {
      // Add the variable if it doesn't already exist
      if (!solution_set.contains(std::make_pair(index, DependencyType::Normal)))
        {
          solution_set[std::make_pair(index, DependencyType::Normal)] =
            std::make_unique<VectorType>();
          solution_transfer_set[std::make_pair(index, DependencyType::Normal)] =
            std::make_unique<SolutionTransfer>(
              matrix_free_handler.get_matrix_free()->get_dof_handler(index));
          new_solution_set[index] = std::make_unique<VectorType>();
        }
      new_solution_set.try_emplace(index, std::make_unique<VectorType>());

      // Add dependencies if they don't exist
      Types::Index field_index = 0;
      for (const auto &dependency_set : variable.get_eval_flag_set_rhs())
        {
          Types::Index dep_index = 0;
          for (const auto &value : dependency_set)
            {
              if (value == dealii::EvaluationFlags::EvaluationFlags::nothing)
                {
                  dep_index++;
                  continue;
                }
              solution_set.try_emplace(
                std::make_pair(field_index, static_cast<DependencyType>(dep_index)),
                std::make_unique<VectorType>());
              solution_transfer_set.try_emplace(
                std::make_pair(field_index, static_cast<DependencyType>(dep_index)),
                std::make_unique<SolutionTransfer>(
                  matrix_free_handler.get_matrix_free()->get_dof_handler(field_index)));

              dep_index++;
            }

          field_index++;
        }
      field_index = 0;
      for (const auto &dependency_set : variable.get_eval_flag_set_lhs())
        {
          Types::Index dep_index = 0;
          for (const auto &value : dependency_set)
            {
              if (value == dealii::EvaluationFlags::EvaluationFlags::nothing)
                {
                  dep_index++;
                  continue;
                }
              solution_set.try_emplace(
                std::make_pair(field_index, static_cast<DependencyType>(dep_index)),
                std::make_unique<VectorType>());
              solution_transfer_set.try_emplace(
                std::make_pair(field_index, static_cast<DependencyType>(dep_index)),
                std::make_unique<SolutionTransfer>(
                  matrix_free_handler.get_matrix_free()->get_dof_handler(field_index)));

              dep_index++;
            }

          field_index++;
        }
    }

  // Initialize the entries according to the corresponding matrix free index
  for (const auto &[pair, solution] : solution_set)
    {
      matrix_free_handler.get_matrix_free()->initialize_dof_vector(*solution, pair.first);
      // TODO (landinjm): Should I ghost values here?
    }
  for (const auto &[index, new_solution] : new_solution_set)
    {
      matrix_free_handler.get_matrix_free()->initialize_dof_vector(*new_solution, index);
      // TODO (landinjm): Should I ghost values here?
    }
}

template <unsigned int dim, typename number>
void
SolutionHandler<dim, number>::reinit(MatrixfreeHandler<dim, number> &matrix_free_handler)
{
  // Initialize the entries according to the corresponding matrix free index
  for (const auto &[pair, solution] : solution_set)
    {
      matrix_free_handler.get_matrix_free()->initialize_dof_vector(*solution, pair.first);
      // TODO (landinjm): Should I ghost values here?
    }
  for (const auto &[index, new_solution] : new_solution_set)
    {
      matrix_free_handler.get_matrix_free()->initialize_dof_vector(*new_solution, index);
      // TODO (landinjm): Should I ghost values here?
    }
}

template <unsigned int dim, typename number>
void
SolutionHandler<dim, number>::mg_init(
  const dealii::MGLevelObject<MatrixfreeHandler<dim, float>> &mg_matrix_free_handler)
{
  // Create all entries and initialize them
  for (unsigned int level = 0; level < mg_solution_set.size(); level++)
    {
      for (unsigned int index = 0; index < mg_solution_set[level].size(); index++)
        {
          mg_solution_set[level][index] = std::make_unique<MGVectorType>();
          mg_matrix_free_handler[level + global_min_level]
            .get_matrix_free()
            ->initialize_dof_vector(*mg_solution_set[level][index], index);
          mg_solution_set[level][index]->update_ghost_values();
        }
    }
}

template <unsigned int dim, typename number>
void
SolutionHandler<dim, number>::mg_reinit(
  const dealii::MGLevelObject<MatrixfreeHandler<dim, float>> &mg_matrix_free_handler)
{
  // Loop over all entries and reinitialize them
  for (unsigned int level = 0; level < mg_solution_set.size(); level++)
    {
      for (unsigned int index = 0; index < mg_solution_set[level].size(); index++)
        {
          mg_matrix_free_handler[level + global_min_level]
            .get_matrix_free()
            ->initialize_dof_vector(*mg_solution_set[level][index], index);
          mg_solution_set[level][index]->update_ghost_values();
        }
    }
}

template <unsigned int dim, typename number>
void
SolutionHandler<dim, number>::update_ghosts() const
{
  for (const auto &[pair, solution] : solution_set)
    {
      solution->update_ghost_values();
    }
  for (const auto &index_vector : mg_solution_set)
    {
      for (const auto &solution : index_vector)
        {
          solution->update_ghost_values();
        }
    }
}

template <unsigned int dim, typename number>
void
SolutionHandler<dim, number>::zero_out_ghosts() const
{
  for (const auto &[index, solution] : new_solution_set)
    {
      solution->zero_out_ghost_values();
    }
}

template <unsigned int dim, typename number>
void
SolutionHandler<dim, number>::apply_constraints(
  unsigned int                             index,
  const dealii::AffineConstraints<number> &constraints)
{
  for (auto &[pair, vector] : solution_set)
    {
      if (pair.first != index)
        {
          continue;
        }
      constraints.distribute(*vector);
    }
}

template <unsigned int dim, typename number>
void
SolutionHandler<dim, number>::apply_initial_condition_for_old_fields()
{
  for (auto &[pair, vector] : solution_set)
    {
      if (pair.second == DependencyType::Normal)
        {
          continue;
        }
      *(get_solution_vector(pair.first, pair.second)) =
        *(get_solution_vector(pair.first, DependencyType::Normal));
    }
}

template <unsigned int dim, typename number>
void
SolutionHandler<dim, number>::update(FieldSolveType field_solve_type,
                                     Types::Index   solve_block,
                                     Types::Index   variable_index)
{
  // Helper function to swap vectors for all dependency types
  auto swap_all_dependency_vectors = [this](Types::Index index, auto &new_vector)
  {
    // Always swap the Normal dependency
    new_vector->swap(*(solution_set.at(std::make_pair(index, DependencyType::Normal))));

    // Swap old dependency types if they exist
    const std::array<DependencyType, 4> old_types = {
      {DependencyType::OldOne,
       DependencyType::OldTwo,
       DependencyType::OldThree,
       DependencyType::OldFour}
    };

    for (const auto &dep_type : old_types)
      {
        if (solution_set.contains(std::make_pair(index, dep_type)))
          {
            new_vector->swap(*(solution_set.at(std::make_pair(index, dep_type))));
          }
      }
  };

  // Loop through the solutions and swap them
  for (auto &[index, new_vector] : new_solution_set)
    {
      const auto &attr_field_type = attributes_list->at(index).get_field_solve_type();

      // Skip if the solve block is wrong
      if (attributes_list->at(index).get_solve_block() != solve_block)
        {
          continue;
        }

      // Skip if field solve types don't match
      if (attr_field_type != field_solve_type)
        {
          continue;
        }

      switch (field_solve_type)
        {
          case FieldSolveType::ExplicitConstant:
            break;
          case FieldSolveType::Explicit:
          case FieldSolveType::ExplicitPostprocess:
            // For ExplicitPostprocess we only swap Normal, but the helper function will
            // do that first and ignore the rest since they should exist
            swap_all_dependency_vectors(index, new_vector);
            break;
          case FieldSolveType::NonexplicitLinear:
            if (variable_index == index)
              {
                swap_all_dependency_vectors(index, new_vector);
                // Additional swap for NonexplicitLinear since the change term is the NEW
                // vector and the Normal vector is the old one
                new_vector->swap(
                  *(solution_set.at(std::make_pair(index, DependencyType::Normal))));
              }
            break;
          case FieldSolveType::NonexplicitAuxiliary:
            if (variable_index == index)
              {
                swap_all_dependency_vectors(index, new_vector);
              }
            break;
          case FieldSolveType::NonexplicitSelfnonlinear:
          case FieldSolveType::NonexplicitCononlinear:
            if (variable_index == index)
              {
                swap_all_dependency_vectors(index, new_vector);

                if (attributes_list->at(index).get_pde_type() != PDEType::Auxiliary)
                  {
                    // Additional swap for NONEXPLICIT_LINEAR since the change term is the
                    // NEW vector and the NORMAL vector is the old one
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, DependencyType::Normal))));
                  }
              }
            break;
          default:
            AssertThrow(false, dealii::ExcMessage("Invalid FieldSolveType"));
        }
    }
}

#include "core/solution_handler.inst"

PRISMS_PF_END_NAMESPACE
