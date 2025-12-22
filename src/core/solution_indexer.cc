// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/solution_indexer.h>

#include <prismspf/config.h>

#include "prismspf/core/group_solution_handler.h"

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
SolutionIndexer<dim, number>::SolutionIndexer(
  unsigned int                                       num_fields, // attributes_list.size()
  const std::set<GroupSolutionHandler<dim, number>> &solution_handlers)
{
  solutions.clear();
  solutions.resize(num_fields, nullptr);
  for (const GroupSolutionHandler<dim, number> &solution_handler : solution_handlers)
    {
      for (unsigned int field_index : solution_handler.get_solve_group().field_indices)
        {
          solutions[field_index] = &solution_handler;
        }
    }
}

template <unsigned int dim, typename number>
auto
SolutionIndexer<dim, number>::get_solution_vector(unsigned int global_index,
                                                  unsigned int relative_level)
  -> const SolutionVector &
{
  return solutions[global_index]->get_solution_vector(global_index, relative_level);
}

template <unsigned int dim, typename number>
auto
SolutionIndexer<dim, number>::get_old_solution_vector(unsigned int age,
                                                      unsigned int global_index,
                                                      unsigned int relative_level)
  -> const SolutionVector &
{
  return solutions[global_index]->get_old_solution_vector(age,
                                                          global_index,
                                                          relative_level);
}

template <unsigned int dim, typename number>
auto
SolutionIndexer<dim, number>::get_new_solution_vector(unsigned int global_index,
                                                      unsigned int relative_level)
  -> const SolutionVector &
{
  return solutions[global_index]->get_new_solution_vector(global_index, relative_level);
}

template <unsigned int dim, typename number>
auto
SolutionIndexer<dim, number>::get_matrix_free(unsigned int global_index,
                                              unsigned int relative_level)
  -> const MatrixFree &
{
  return solutions[global_index]->get_matrix_free(relative_level);
}

template <unsigned int dim, typename number>
auto
SolutionIndexer<dim, number>::get_solution_level_and_block_index(
  unsigned int global_index,
  unsigned int relative_level)
  -> std::pair<const SolutionLevel<dim, number> &, unsigned int>
{
  auto &sol = solutions[global_index];
  return {sol->get_solution_level(relative_level), sol->get_block_index(global_index)};
}

template <unsigned int dim, typename number>
unsigned int
SolutionIndexer<dim, number>::get_block_index(unsigned int global_index) const

{
  return solutions[global_index]->get_block_index(global_index);
}

// #include "core/solution_indexer.inst"

PRISMS_PF_END_NAMESPACE
