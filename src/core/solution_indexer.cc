// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/solution_indexer.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

// TODO add initialization checks.

template <unsigned int dim, typename number>
void
SolutionIndexer<dim, number>::init(
  unsigned int                                            num_fields,
  const std::vector<GroupSolutionHandler<dim, number> *> &solution_handlers)
{
  solutions.clear();
  solutions.resize(num_fields, nullptr);
  block_indices.clear();
  block_indices.resize(num_fields, -1);
  for (GroupSolutionHandler<dim, number> *solution_handler : solution_handlers)
    {
      for (unsigned int field_index : solution_handler->get_solve_block().field_indices)
        {
          solutions[field_index]     = solution_handler;
          block_indices[field_index] = solution_handler->get_block_index(field_index);
        }
    }
}

template <unsigned int dim, typename number>
auto
SolutionIndexer<dim, number>::get_solution_vector(unsigned int global_index) const
  -> const SolutionVector<number> &
{
  return solutions[global_index]->get_solution_vector(global_index);
}

template <unsigned int dim, typename number>
auto
SolutionIndexer<dim, number>::get_solution_vector(unsigned int global_index)
  -> SolutionVector<number> &
{
  return solutions[global_index]->get_solution_vector(global_index);
}

template <unsigned int dim, typename number>
auto
SolutionIndexer<dim, number>::get_old_solution_vector(unsigned int age,
                                                      unsigned int global_index) const
  -> const SolutionVector<number> &
{
  return solutions[global_index]->get_old_solution_vector(age, global_index);
}

template <unsigned int dim, typename number>
auto
SolutionIndexer<dim, number>::get_old_solution_vector(unsigned int age,
                                                      unsigned int global_index)
  -> SolutionVector<number> &
{
  return solutions[global_index]->get_old_solution_vector(age, global_index);
}

template <unsigned int dim, typename number>
const SolveBlock &
SolutionIndexer<dim, number>::get_solve_block(unsigned int index) const
{
  return solutions[index]->get_solve_block();
}

template <unsigned int dim, typename number>
auto
SolutionIndexer<dim, number>::get_solution_level_and_block_index(
  unsigned int global_index,
  unsigned int relative_level) const
  -> std::pair<const SolutionLevel<dim, number> *, unsigned int>
{
  auto &sol = solutions[global_index];
  if (relative_level == -1)
    {
      return {&(sol->get_primary_solutions()), sol->get_block_index(global_index)};
    }
  // else
  return {&(sol->get_solution_level(relative_level)), sol->get_block_index(global_index)};
}

template <unsigned int dim, typename number>
unsigned int
SolutionIndexer<dim, number>::get_block_index(unsigned int global_index) const
{
  // return solutions[global_index]->get_block_index(global_index);
  return block_indices[global_index];
}

template <unsigned int dim, typename number>
std::vector<unsigned int>
SolutionIndexer<dim, number>::get_block_indices() const
{
  return block_indices;
}

#include "core/solution_indexer.inst"

PRISMS_PF_END_NAMESPACE
