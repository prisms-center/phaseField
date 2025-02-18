#include <deal.II/lac/la_parallel_vector.h>

#include <prismspf/config.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <unordered_map>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
solutionHandler<dim>::solutionHandler(
  const std::map<unsigned int, variableAttributes> &_attributes_list)
  : attributes_list(_attributes_list)
{}

template <int dim>
solutionHandler<dim>::~solutionHandler()
{
  for (const auto &[pair, solution] : solution_set)
    {
      delete solution;
    }
  solution_set.clear();
}

template <int dim>
void
solutionHandler<dim>::init(matrixfreeHandler<dim> &matrix_free_handler)
{
  // Create all entries
  for (const auto &[index, variable] : attributes_list)
    {
      // Add the current variable if it doesn't already exist
      if (solution_set.find(std::make_pair(index, dependencyType::NORMAL)) ==
          solution_set.end())
        {
          solution_set[std::make_pair(index, dependencyType::NORMAL)] = new VectorType();

          new_solution_set[index] = new VectorType();
        }
      else if (new_solution_set.find(index) == new_solution_set.end())
        {
          new_solution_set[index] = new VectorType();
        }

      // Add dependencies if they don't exist
      for (const auto &[pair, flags] : variable.eval_flag_set_RHS)
        {
          if (solution_set.find(pair) == solution_set.end())
            {
              solution_set[pair] = new VectorType();
            }
        }
      for (const auto &[pair, flags] : variable.eval_flag_set_LHS)
        {
          if (solution_set.find(pair) == solution_set.end())
            {
              solution_set[pair] = new VectorType();
            }
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
            if (attributes_list.at(index).field_solve_type == field_solve_type)
              {
                (*new_vector)
                  .swap(
                    *(solution_set.at(std::make_pair(index, dependencyType::NORMAL))));
                if (solution_set.find(std::make_pair(index, dependencyType::OLD_1)) !=
                    solution_set.end())
                  {
                    (*new_vector)
                      .swap(
                        *(solution_set.at(std::make_pair(index, dependencyType::OLD_1))));
                  }
                if (solution_set.find(std::make_pair(index, dependencyType::OLD_2)) !=
                    solution_set.end())
                  {
                    (*new_vector)
                      .swap(
                        *(solution_set.at(std::make_pair(index, dependencyType::OLD_2))));
                  }
                if (solution_set.find(std::make_pair(index, dependencyType::OLD_3)) !=
                    solution_set.end())
                  {
                    (*new_vector)
                      .swap(
                        *(solution_set.at(std::make_pair(index, dependencyType::OLD_3))));
                  }
                if (solution_set.find(std::make_pair(index, dependencyType::OLD_4)) !=
                    solution_set.end())
                  {
                    (*new_vector)
                      .swap(
                        *(solution_set.at(std::make_pair(index, dependencyType::OLD_4))));
                  }
              }
            break;
          case fieldSolveType::NONEXPLICIT_LINEAR:
            if (attributes_list.at(index).field_solve_type == field_solve_type &&
                variable_index == index)
              {
                (*new_vector)
                  .swap(
                    *(solution_set.at(std::make_pair(index, dependencyType::NORMAL))));
                if (solution_set.find(std::make_pair(index, dependencyType::OLD_1)) !=
                    solution_set.end())
                  {
                    (*new_vector)
                      .swap(
                        *(solution_set.at(std::make_pair(index, dependencyType::OLD_1))));
                  }
                if (solution_set.find(std::make_pair(index, dependencyType::OLD_2)) !=
                    solution_set.end())
                  {
                    (*new_vector)
                      .swap(
                        *(solution_set.at(std::make_pair(index, dependencyType::OLD_2))));
                  }
                if (solution_set.find(std::make_pair(index, dependencyType::OLD_3)) !=
                    solution_set.end())
                  {
                    (*new_vector)
                      .swap(
                        *(solution_set.at(std::make_pair(index, dependencyType::OLD_3))));
                  }
                if (solution_set.find(std::make_pair(index, dependencyType::OLD_4)) !=
                    solution_set.end())
                  {
                    (*new_vector)
                      .swap(
                        *(solution_set.at(std::make_pair(index, dependencyType::OLD_4))));
                  }
                (*new_vector)
                  .swap(
                    *(solution_set.at(std::make_pair(index, dependencyType::NORMAL))));
              }
            break;
          case fieldSolveType::NONEXPLICIT_SELF_NONLINEAR:
            break;
          case fieldSolveType::NONEXPLICIT_AUXILIARY:
            if (attributes_list.at(index).field_solve_type == field_solve_type &&
                variable_index == index)
              {
                (*new_vector)
                  .swap(
                    *(solution_set.at(std::make_pair(index, dependencyType::NORMAL))));
                if (solution_set.find(std::make_pair(index, dependencyType::OLD_1)) !=
                    solution_set.end())
                  {
                    (*new_vector)
                      .swap(
                        *(solution_set.at(std::make_pair(index, dependencyType::OLD_1))));
                  }
                if (solution_set.find(std::make_pair(index, dependencyType::OLD_2)) !=
                    solution_set.end())
                  {
                    (*new_vector)
                      .swap(
                        *(solution_set.at(std::make_pair(index, dependencyType::OLD_2))));
                  }
                if (solution_set.find(std::make_pair(index, dependencyType::OLD_3)) !=
                    solution_set.end())
                  {
                    (*new_vector)
                      .swap(
                        *(solution_set.at(std::make_pair(index, dependencyType::OLD_3))));
                  }
                if (solution_set.find(std::make_pair(index, dependencyType::OLD_4)) !=
                    solution_set.end())
                  {
                    (*new_vector)
                      .swap(
                        *(solution_set.at(std::make_pair(index, dependencyType::OLD_4))));
                  }
              }
            break;
          case fieldSolveType::NONEXPLICIT_CO_NONLINEAR:
            break;
          case fieldSolveType::EXPLICIT_POSTPROCESS:
            if (attributes_list.at(index).field_solve_type == field_solve_type)
              {
                (*new_vector)
                  .swap(
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