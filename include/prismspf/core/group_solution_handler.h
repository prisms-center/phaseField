// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/constraint_manager.h>
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/matrix_free_manager.h>
#include <prismspf/core/solve_block.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief The solution vectors with respect to on some multigrid level.
 *
 * `solutions` is a block vector of all the fields that we keep track of at this level.
 * `old_solutions` is a vector of with length equal to the number of old states we keep
 * track old. Importantly, for each std::vector entry the block vector has the same shape
 * as `solutions`.
 *
 * @remark There is a [minor] optimization we could make regarding the old solutions. If a
 * user wants to store up the 2nd old state, we allocate up to the 2nd old state for all
 * fields. This behavior is indiscriminate on whether the user actually needs those
 * vectors, so we could not do that. However, it's painful to deal with and will likely
 * only show minor performance improvements.
 */
template <unsigned int dim, typename number>
struct SolutionLevel
{
  BlockVector<number>              solutions;
  std::vector<BlockVector<number>> old_solutions;
};

/**
 * @brief Class that manages solution initialization and swapping with old solutions.
 */
template <unsigned int dim, typename number>
class GroupSolutionHandler
{
public:
#if DEAL_II_VERSION_MAJOR >= 9 && DEAL_II_VERSION_MINOR >= 7
  using SolutionTransfer = dealii::SolutionTransfer<dim, SolutionVector<number>>;
#else
  using SolutionTransfer =
    dealii::parallel::distributed::SolutionTransfer<dim, SolutionVector<number>>;
#endif

  /**
   * @brief Constructor.
   */
  GroupSolutionHandler(SolveBlock                            _solve_block,
                       const std::vector<FieldAttributes>   &_attributes_list,
                       const MatrixFreeManager<dim, number> &_matrix_free_manager);

  /**
   * @brief Get the solution vector set. This contains all the normal fields and is
   * typically used for output.
   */
  [[nodiscard]] BlockVector<number> &
  get_solution_full_vector(unsigned int relative_level = 0);

  /**
   * @brief Get the const solution vector set. This contains all the normal fields and is
   * typically used for output.
   */
  [[nodiscard]] const BlockVector<number> &
  get_solution_full_vector(unsigned int relative_level = 0) const;

  /**
   * @brief Get a solution vector of a given field index.
   */
  [[nodiscard]] SolutionVector<number> &
  get_solution_vector(unsigned int global_index, unsigned int relative_level = 0);

  /**
   * @brief Get a solution vector of a given field index.
   */
  [[nodiscard]] const SolutionVector<number> &
  get_solution_vector(unsigned int global_index, unsigned int relative_level = 0) const;

  /**
   * @brief Get the old solution vector set at a given age.
   */
  [[nodiscard]] BlockVector<number> &
  get_old_solution_full_vector(unsigned int age, unsigned int relative_level = 0);

  /**
   * @brief Get the old solution vector set at a given age.
   */
  [[nodiscard]] const BlockVector<number> &
  get_old_solution_full_vector(unsigned int age, unsigned int relative_level = 0) const;

  /**
   * @brief Get a solution vector of a given field index at a given age.
   */
  [[nodiscard]] SolutionVector<number> &
  get_old_solution_vector(unsigned int age,
                          unsigned int global_index,
                          unsigned int relative_level = 0);

  /**
   * @brief Get a solution vector of a given field index at a given age.
   */
  [[nodiscard]] const SolutionVector<number> &
  get_old_solution_vector(unsigned int age,
                          unsigned int global_index,
                          unsigned int relative_level = 0) const;

  /**
   * @brief Get the solutions object at a level.
   */
  [[nodiscard]] SolutionLevel<dim, number> &
  get_solution_level(unsigned int relative_level = 0);

  /**
   * @brief Get the solutions object at a level.
   */
  [[nodiscard]] const SolutionLevel<dim, number> &
  get_solution_level(unsigned int relative_level = 0) const;

  /**
   * @brief Get the underlying matrix_free objects.
   */
  const std::vector<MatrixFree<dim, number>> &
  get_matrix_free_levels() const;

  /**
   * @brief Get the block index from the global field index.
   */
  [[nodiscard]] unsigned int
  get_block_index(unsigned int global_index) const;

  /**
   * @brief Get the underlying solve block object.
   */
  [[nodiscard]] const SolveBlock &
  get_solve_block() const;

  /**
   * @brief Get the block index from the global index.
   */
  [[nodiscard]] const std::vector<unsigned int> &
  get_global_to_block_index() const;

  /**
   * @brief Get the global index from the block index.
   */
  [[nodiscard]] const std::vector<unsigned int> &
  get_block_to_global_index() const;

  /**
   * @brief Initialize the solution set.
   * @pre MatrixFree is already reinit.
   */
  void
  init(const std::vector<unsigned int> &max_age_per_level);

  /**
   * @brief Reinitialize the solution set.
   * @pre MatrixFree is already reinit.
   */
  void
  reinit();

  /**
   * @brief Update the ghost values.
   */
  void
  update_ghosts(unsigned int relative_level = 0) const;

  /**
   * @brief Zero out the ghost values.
   */
  void
  zero_out_ghosts(unsigned int relative_level = 0) const;

  /**
   * @brief Apply the constraints to the solution vector.
   */
  void
  apply_constraints(unsigned int relative_level = 0);

  /**
   * @brief Apply the constraints to another vector.
   */
  void
  apply_constraints(BlockVector<number> &solution_vector,
                    unsigned int         relative_level = 0);

  /**
   * @brief Apply the given constraints to all the solution vectors, including old.
   */
  void
  apply_constraints_to_all(unsigned int relative_level);

  /**
   * @brief Apply initial condition to the old fields. For now, this simply copies the
   * values in the normal field to the old.
   */
  void
  apply_initial_condition_for_old_fields();

  /**
   * @brief Update and propagate the old solutions.
   */
  void
  update(unsigned int relative_level = 0);

  /**
   * @brief Reinit the solution transfer objects.
   */
  void
  init_solution_transfer();

  /**
   * @brief Prepare for solution transfer
   */
  void
  prepare_for_solution_transfer();

  /**
   * @brief Transfer solutions
   */
  void
  execute_solution_transfer();

  /**
   * @brief Reinit the solution transfer objects.
   */
  template <unsigned int degree>
  void
  reinit_mg_transfer(const DoFManager<dim, degree> &dof_manager);

  /**
   * @brief Transfer solutions to mg levels
   */
  template <unsigned int degree>
  void
  mg_transfer_down(const DoFManager<dim, degree> &dof_manager,
                   unsigned int                   finest_level,
                   bool                           transfer_old_solutions = false);

  /**
   * @brief Number of refinement levels that will be tracked.
   */
  [[nodiscard]] unsigned int
  num_levels();

  /**
   * @brief Print the solution vector set.
   */
  void
  print_solution_full_vector(std::ostream &out, unsigned int relative_level = 0) const;

private:
  /**
   * @brief Information about the solve block this handler is responsible for.
   */
  SolveBlock solve_block;

  /**
   * @brief Mapping from block index to global field index.
   */
  std::vector<unsigned int> block_to_global_index;

  /**
   * @brief Mapping from global field index to block index.
   */
  std::vector<unsigned int> global_to_block_index;

  /**
   * @brief Solutions of each level.
   */
  std::vector<SolutionLevel<dim, number>> solution_levels;

  /**
   * @brief Pointer to MatrixFree manager
   */
  const MatrixFreeManager<dim, number> *matrix_free_manager = nullptr;

  /**
   * @brief Utility for solution transfer to different mesh (for AMR). Can only work on
   * one block at a time.
   * @note solution transfers can work on multiple solutions as long as they are using
   * the same underlying dof numbering. All of our scalars and vectors have identical
   * dof handlers, so we may be able to take advantage of this instead of creating
   * several solution transfers as we do here. I don't know if this affects performance.
   * todo
   */
  std::vector<SolutionTransfer> block_solution_transfer;

  /**
   * @brief Constraints needed for mg transfer object.
   */
  std::vector<dealii::MGConstrainedDoFs> mg_constraints;

  /**
   * @brief Utility object to transfer solutions between multigrid levels.
   */
  std::vector<dealii::MGTransferMatrixFree<dim, number>> mg_transfer;
};

template <unsigned int dim, typename number>
template <unsigned int degree>
inline void
GroupSolutionHandler<dim, number>::reinit_mg_transfer(
  const DoFManager<dim, degree> &dof_manager)
{
  // 1. Initialize constraints
  const unsigned int num_blocks = solve_block.field_indices.size();
  mg_constraints                = std::vector<dealii::MGConstrainedDoFs>(num_blocks);
  mg_transfer = std::vector<dealii::MGTransferMatrixFree<dim, number>>(num_blocks);
  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
    {
      unsigned int field_index = block_to_global_index[block_index];
      mg_constraints[block_index].initialize(
        dof_manager.get_field_dof_handler(field_index, 0));
      ConditionalOStreams::pout_base()
        << dof_manager.get_field_dof_handler(field_index, 0).n_dofs();
      // 2. Initialize MG Transfer
      mg_transfer[block_index].initialize_constraints(mg_constraints[block_index]);
      mg_transfer[block_index].build(dof_manager.get_field_dof_handler(field_index, 0));
    }
}

template <unsigned int dim, typename number>
template <unsigned int degree>
inline void
GroupSolutionHandler<dim, number>::mg_transfer_down(
  const DoFManager<dim, degree> &dof_manager,
  unsigned int                   finest_level,
  bool                           transfer_old_solutions)
{
  const unsigned int num_blocks = solve_block.field_indices.size();
  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
    {
      unsigned int field_index = block_to_global_index[block_index];
      dealii::MGLevelObject<SolutionVector<number>> temp_mg_solutions(
        1 + finest_level - solution_levels.size(),
        finest_level);

      const auto &dof_handler = dof_manager.get_field_dof_handler(field_index, 0);

      // transfer regular solutions to mg levels
      mg_transfer[block_index].interpolate_to_mg(dof_handler,
                                                 temp_mg_solutions,
                                                 solution_levels[0].solutions.block(
                                                   block_index));
      // swap to actual mg solution vectors
      for (unsigned int relative_level = 0; relative_level < solution_levels.size();
           ++relative_level)
        {
          unsigned int level = finest_level - relative_level;
          ConditionalOStreams::pout_base()
            << "R: " << relative_level << "\nL: " << level << "\n Pre: "
            << solution_levels[relative_level].solutions.block(block_index).size();
          solution_levels[relative_level]
            .solutions.block(block_index)
            .swap(temp_mg_solutions[level]);
          ConditionalOStreams::pout_base()
            << "\n Post: "
            << solution_levels[relative_level].solutions.block(block_index).size() << "\n"
            << std::flush;
        }

      if (!transfer_old_solutions)
        {
          return;
        }
      // transfer old solutions
      for (unsigned int age_index = 0;
           age_index < solution_levels[0].old_solutions.size();
           ++age_index)
        {
          mg_transfer[block_index].interpolate_to_mg(
            dof_handler,
            temp_mg_solutions,
            solution_levels[0].old_solutions[age_index].block(block_index));
          for (unsigned int relative_level = 0; relative_level < solution_levels.size();
               ++relative_level)
            {
              unsigned int level = finest_level - relative_level;
              if (age_index < solution_levels[relative_level].old_solutions.size())
                {
                  solution_levels[relative_level]
                    .old_solutions[age_index]
                    .block(block_index)
                    .swap(temp_mg_solutions[level]);
                }
            }
        }
    }
}

PRISMS_PF_END_NAMESPACE
