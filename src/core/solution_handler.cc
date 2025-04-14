// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/lac/affine_constraints.h>

#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <map>
#include <memory>
#include <string>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
solutionHandler<dim>::solutionHandler(
  const std::map<unsigned int, variableAttributes> &_attributes_list)
  : attributes_list(&_attributes_list)
{}

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
void
solutionHandler<dim>::init(matrixfreeHandler<dim> &matrix_free_handler)
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
    }
  for (const auto &[index, new_solution] : new_solution_set)
    {
      matrix_free_handler.get_matrix_free()->initialize_dof_vector(*new_solution, index);
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
  // Loop through the solutions and swap them
  for (auto &[index, new_vector] : new_solution_set)
    {
      switch (field_solve_type)
        {
          case fieldSolveType::EXPLICIT_CONSTANT:
            break;
          case fieldSolveType::EXPLICIT:
            if (attributes_list->at(index).field_solve_type == field_solve_type)
              {
                new_vector->swap(
                  *(solution_set.at(std::make_pair(index, dependencyType::NORMAL))));
                if (solution_set.contains(std::make_pair(index, dependencyType::OLD_1)))
                  {
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, dependencyType::OLD_1))));
                  }
                if (solution_set.contains(std::make_pair(index, dependencyType::OLD_2)))
                  {
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, dependencyType::OLD_2))));
                  }
                if (solution_set.contains(std::make_pair(index, dependencyType::OLD_3)))
                  {
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, dependencyType::OLD_3))));
                  }
                if (solution_set.contains(std::make_pair(index, dependencyType::OLD_4)))
                  {
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, dependencyType::OLD_4))));
                  }
              }
            break;
          case fieldSolveType::NONEXPLICIT_LINEAR:
            if (attributes_list->at(index).field_solve_type == field_solve_type &&
                variable_index == index)
              {
                new_vector->swap(
                  *(solution_set.at(std::make_pair(index, dependencyType::NORMAL))));
                if (solution_set.contains(std::make_pair(index, dependencyType::OLD_1)))
                  {
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, dependencyType::OLD_1))));
                  }
                if (solution_set.contains(std::make_pair(index, dependencyType::OLD_2)))
                  {
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, dependencyType::OLD_2))));
                  }
                if (solution_set.contains(std::make_pair(index, dependencyType::OLD_3)))
                  {
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, dependencyType::OLD_3))));
                  }
                if (solution_set.contains(std::make_pair(index, dependencyType::OLD_4)))
                  {
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, dependencyType::OLD_4))));
                  }
                new_vector->swap(
                  *(solution_set.at(std::make_pair(index, dependencyType::NORMAL))));
              }
            break;
          case fieldSolveType::NONEXPLICIT_SELF_NONLINEAR:
            break;
          case fieldSolveType::NONEXPLICIT_AUXILIARY:
            if (attributes_list->at(index).field_solve_type == field_solve_type &&
                variable_index == index)
              {
                new_vector->swap(
                  *(solution_set.at(std::make_pair(index, dependencyType::NORMAL))));
                if (solution_set.contains(std::make_pair(index, dependencyType::OLD_1)))
                  {
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, dependencyType::OLD_1))));
                  }
                if (solution_set.contains(std::make_pair(index, dependencyType::OLD_2)))
                  {
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, dependencyType::OLD_2))));
                  }
                if (solution_set.contains(std::make_pair(index, dependencyType::OLD_3)))
                  {
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, dependencyType::OLD_3))));
                  }
                if (solution_set.contains(std::make_pair(index, dependencyType::OLD_4)))
                  {
                    new_vector->swap(
                      *(solution_set.at(std::make_pair(index, dependencyType::OLD_4))));
                  }
              }
            break;
          case fieldSolveType::NONEXPLICIT_CO_NONLINEAR:
            break;
          case fieldSolveType::EXPLICIT_POSTPROCESS:
            if (attributes_list->at(index).field_solve_type == field_solve_type)
              {
                new_vector->swap(
                  *(solution_set.at(std::make_pair(index, dependencyType::NORMAL))));
              }
            break;
          default:
            AssertThrow(false, dealii::ExcMessage("Invalid fieldSolveType"));
        }
    }
}

INSTANTIATE_UNI_TEMPLATE(solutionHandler)

PRISMS_PF_END_NAMESPACE
