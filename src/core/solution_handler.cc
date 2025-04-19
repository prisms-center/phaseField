// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/lac/affine_constraints.h>

#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <array>
#include <map>
#include <memory>
#include <string>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
solutionHandler<dim>::solutionHandler(
  const std::map<unsigned int, variableAttributes> &_attributes_list,
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

template <int dim>
std::map<unsigned int, typename solutionHandler<dim>::VectorType *>
solutionHandler<dim>::get_solution_vector() const
{
  std::map<unsigned int, VectorType *> temp;
  for (const auto &[pair, vector] : solution_set)
    {
      if (pair.second != dependencyType::NORMAL)
        {
          continue;
        }

      Assert(vector != nullptr, dealii::ExcNotInitialized());

      temp.emplace(pair.first, vector.get());
    }

  return temp;
}

template <int dim>
typename solutionHandler<dim>::VectorType *
solutionHandler<dim>::get_solution_vector(unsigned int   index,
                                          dependencyType dependency_type) const
{
  const auto pair = std::make_pair(index, dependency_type);

  Assert(solution_set.contains(pair),
         dealii::ExcMessage(
           "There is no solution vector for the given index = " + std::to_string(index) +
           " and type = " + to_string(dependency_type)));
  Assert(solution_set.at(pair) != nullptr, dealii::ExcNotInitialized());

  return solution_set.at(pair).get();
}

template <int dim>
std::map<unsigned int, typename solutionHandler<dim>::VectorType *>
solutionHandler<dim>::get_new_solution_vector() const
{
  std::map<unsigned int, VectorType *> temp;
  for (const auto &[index, vector] : new_solution_set)
    {
      Assert(vector != nullptr, dealii::ExcNotInitialized());

      temp.emplace(index, vector.get());
    }

  return temp;
}

template <int dim>
typename solutionHandler<dim>::VectorType *
solutionHandler<dim>::get_new_solution_vector(unsigned int index) const
{
  Assert(new_solution_set.contains(index),
         dealii::ExcMessage("There is no new solution vector for the given index = " +
                            std::to_string(index)));
  Assert(new_solution_set.at(index) != nullptr, dealii::ExcNotInitialized());

  return new_solution_set.at(index).get();
}

template <int dim>
std::vector<typename solutionHandler<dim>::MGVectorType *>
solutionHandler<dim>::get_mg_solution_vector(unsigned int level) const
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

template <int dim>
typename solutionHandler<dim>::MGVectorType *
solutionHandler<dim>::get_mg_solution_vector(unsigned int level, unsigned int index) const
{
  Assert(has_multigrid, dealii::ExcNotInitialized());

  // Convert absolute level to a relative level
  const unsigned int relative_level = level - global_min_level;
  Assert(relative_level < mg_solution_set.size(),
         dealii::ExcMessage("The mg solution set does not contain level = " +
                            std::to_string(level)));
  const unsigned int global_index = mg_info->get_global_index(index, relative_level);
  Assert(index < mg_solution_set[relative_level].size(),
         dealii::ExcMessage(
           "The mg solution at the given level does not contain index = " +
           std::to_string(global_index)));
  return mg_solution_set[relative_level][index].get();
}

template <int dim>
void
solutionHandler<dim>::init(matrixfreeHandler<dim, double> &matrix_free_handler)
{
  // Create all entries
  for (const auto &[index, variable] : *attributes_list)
    {
      // Add the current variable if it doesn't already exist
      if (!solution_set.contains(std::make_pair(index, dependencyType::NORMAL)))
        {
          solution_set[std::make_pair(index, dependencyType::NORMAL)] =
            std::make_unique<VectorType>();

          new_solution_set[index] = std::make_unique<VectorType>();
        }
      new_solution_set.try_emplace(index, std::make_unique<VectorType>());

      // Add dependencies if they don't exist
      for (const auto &[pair, flags] : variable.eval_flag_set_RHS)
        {
          solution_set.try_emplace(pair, std::make_unique<VectorType>());
        }
      for (const auto &[pair, flags] : variable.eval_flag_set_LHS)
        {
          solution_set.try_emplace(pair, std::make_unique<VectorType>());
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

template <int dim>
void
solutionHandler<dim>::mg_init(
  const dealii::MGLevelObject<matrixfreeHandler<dim, float>> &mg_matrix_free_handler)
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

template <int dim>
void
solutionHandler<dim>::update_ghosts() const
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

template <int dim>
void
solutionHandler<dim>::apply_constraints(
  unsigned int                             index,
  const dealii::AffineConstraints<double> &constraints)
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

template <int dim>
void
solutionHandler<dim>::apply_initial_condition_for_old_fields()
{
  for (auto &[pair, vector] : solution_set)
    {
      if (pair.second == dependencyType::NORMAL)
        {
          continue;
        }
      *(get_solution_vector(pair.first, pair.second)) =
        *(get_solution_vector(pair.first, dependencyType::NORMAL));
    }
}

template <int dim>
void
solutionHandler<dim>::update(const fieldSolveType &field_solve_type,
                             const unsigned int   &variable_index)
{
  // Helper function to swap vectors for all dependency types
  auto swap_all_dependency_vectors = [this](unsigned int index, auto &new_vector)
  {
    // Always swap the NORMAL dependency
    new_vector->swap(*(solution_set.at(std::make_pair(index, dependencyType::NORMAL))));

    // Swap old dependency types if they exist
    const std::array<dependencyType, 4> old_types = {dependencyType::OLD_1,
                                                     dependencyType::OLD_2,
                                                     dependencyType::OLD_3,
                                                     dependencyType::OLD_4};

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
      const auto &attr_field_type = attributes_list->at(index).field_solve_type;

      // Skip if field solve types don't match
      if (attr_field_type != field_solve_type)
        {
          continue;
        }

      switch (field_solve_type)
        {
          case fieldSolveType::EXPLICIT_CONSTANT:
            break;
          case fieldSolveType::EXPLICIT:
          case fieldSolveType::EXPLICIT_POSTPROCESS:
            // For EXPLICIT_POSTPROCESS we only swap NORMAL, but the helper function will
            // do that first and ignore the rest since they should exist
            swap_all_dependency_vectors(index, new_vector);
            break;
          case fieldSolveType::NONEXPLICIT_LINEAR:
            if (variable_index == index)
              {
                swap_all_dependency_vectors(index, new_vector);
                // Additional swap for NONEXPLICIT_LINEAR since the change term is the NEW
                // vector and the NORMAL vector is the old one
                new_vector->swap(
                  *(solution_set.at(std::make_pair(index, dependencyType::NORMAL))));
              }
            break;
          case fieldSolveType::NONEXPLICIT_AUXILIARY:
            if (variable_index == index)
              {
                swap_all_dependency_vectors(index, new_vector);
              }
            break;
          case fieldSolveType::NONEXPLICIT_SELF_NONLINEAR:
          case fieldSolveType::NONEXPLICIT_CO_NONLINEAR:
            Assert(false, dealii::ExcNotImplemented());
            break;
          default:
            AssertThrow(false, dealii::ExcMessage("Invalid fieldSolveType"));
        }
    }
}

INSTANTIATE_UNI_TEMPLATE(solutionHandler)

PRISMS_PF_END_NAMESPACE
