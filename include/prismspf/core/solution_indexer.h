// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/group_solution_handler.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Class that provides access to solution vectors spread across different
 * groups.
 */
template <unsigned int dim, typename number>
class SolutionIndexer
{
public:
  /**
   * @brief Constructor.
   */
  SolutionIndexer(unsigned int num_fields, // attributes_list.size()
                  std::vector<GroupSolutionHandler<dim, number> *> solution_handlers);

  /**
   * @brief Get a solution vector of a given field index.
   */
  [[nodiscard]] const SolutionVector<number> &
  get_solution_vector(unsigned int global_index, unsigned int relative_level = 0) const;
  /**
   * @brief Get a solution vector of a given field index.
   */
  [[nodiscard]] SolutionVector<number> &
  get_solution_vector(unsigned int global_index, unsigned int relative_level = 0);

  /**
   * @brief Get a solution vector of a given field index at a given age.
   */
  [[nodiscard]] const SolutionVector<number> &
  get_old_solution_vector(unsigned int age,
                          unsigned int global_index,
                          unsigned int relative_level = 0) const;
  /**
   * @brief Get a solution vector of a given field index at a given age.
   */
  [[nodiscard]] SolutionVector<number> &
  get_old_solution_vector(unsigned int age,
                          unsigned int global_index,
                          unsigned int relative_level = 0);

  /**
   * @brief Get the matrixfree object of the group a given field index.
   */
  [[nodiscard]] const MatrixFree<dim, number> &
  get_matrix_free(unsigned int index, unsigned int relative_level = 0) const;
  /**
   * @brief Get the matrixfree object of the group a given field index.
   */
  [[nodiscard]] MatrixFree<dim, number> &
  get_matrix_free(unsigned int index, unsigned int relative_level = 0);

  /**
   * @brief Get the solve group of a given field index.
   */
  [[nodiscard]] const SolveBlock &
  get_solve_block(unsigned int index) const;

  /**
   * @brief Get the matrixfree object of the group a given field index.
   */
  [[nodiscard]] std::pair<const SolutionLevel<dim, number> *, unsigned int>
  get_solution_level_and_block_index(unsigned int index,
                                     unsigned int relative_level = 0) const;

  /**
   * @brief Get the block index of a field within its solve group
   */
  [[nodiscard]] unsigned int
  get_block_index(unsigned int global_index) const;

private:
  std::vector<GroupSolutionHandler<dim, number> *> solutions;
};

PRISMS_PF_END_NAMESPACE
