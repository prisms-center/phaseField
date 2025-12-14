// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
class MatrixFreeContainer;

template <unsigned int dim>
class MGInfo;

/**
 * @brief Class that manages solution initialization and swapping with old solutions.
 */
template <unsigned int dim, typename number>
class GroupSolutionHandler
{
public:
  using BlockVector      = dealii::LinearAlgebra::distributed::BlockVector<number>;
  using SolutionVector   = BlockVector::BlockType;
  using MGBlockVector    = dealii::LinearAlgebra::distributed::BlockVector<float>;
  using MGSolutionVector = MGBlockVector::BlockType;
#if DEAL_II_VERSION_MAJOR >= 9 && DEAL_II_VERSION_MINOR >= 7
  using SolutionTransfer = dealii::SolutionTransfer<dim, BlockVector>;
#else
  using SolutionTransfer =
    dealii::parallel::distributed::SolutionTransfer<dim, BlockVector>;
#endif

  /**
   * @brief Constructor.
   */
  GroupSolutionHandler(const SolveGroup                   &_solve_group,
                       const std::vector<FieldAttributes> &_attributes_list,
                       const MGInfo<dim>                  &_mg_info);

  // TODO (fractalsbyx): add 'const T& get() const' versions.
  /**
   * @brief Get the solution vector set. This contains all the normal fields and is
   * typically used for output.
   */
  [[nodiscard]] BlockVector &
  get_solution_full_vector(unsigned int relative_level = 0);

  /**
   * @brief Get a solution vector of a given field index.
   */
  [[nodiscard]] SolutionVector &
  get_solution_vector(unsigned int global_index, unsigned int relative_level = 0);

  /**
   * @brief Get the old solution vector set at a given age.
   */
  [[nodiscard]] BlockVector &
  get_old_solution_full_vector(unsigned int age, unsigned int relative_level = 0);

  /**
   * @brief Get a solution vector of a given field index at a given age.
   */
  [[nodiscard]] SolutionVector &
  get_old_solution_vector(unsigned int age,
                          unsigned int global_index,
                          unsigned int relative_level = 0);

  /**
   * @brief Get the "new" solution vector set.
   */
  [[nodiscard]] BlockVector &
  get_new_solution_full_vector(unsigned int relative_level = 0);

  /**
   * @brief Get the "new" solution vector of a given field index.
   */
  [[nodiscard]] SolutionVector &
  get_new_solution_vector(unsigned int index, unsigned int relative_level = 0);

  /**
   * @brief Initialize the solution set.
   */
  void
  init();

  /**
   * @brief Reinitialize the solution set.
   */
  void
  reinit();

  /**
   * @brief Update the ghost values.
   *
   * TODO (landinjm): Fix so this isn't as wasteful in updating ghost values for all
   * solution vectors.
   */
  void
  update_ghosts() const;

  /**
   * @brief Zero out the ghost values.
   *
   * TODO (landinjm): Fix so this isn't as wasteful in zeroing ghost values for all
   * solution vectors.
   */
  void
  zero_out_ghosts() const;

  /**
   * @brief Apply the given constraints to a solution vector of a given field index.
   *
   * Note this applies constraints for all dependencyTypes of the given index.
   */
  void
  apply_constraints(unsigned int                             field_index,
                    const dealii::AffineConstraints<number> &constraints);

  /**
   * @brief Apply intial condition to the old fields. For now, this simply copies the
   * values in the normal field to the old.
   *
   * TODO (landinjm): What should we do for the initial condition of old fields.
   */
  void
  apply_initial_condition_for_old_fields();

  /**
   * @brief Update the `solution_set` with the `new_solution_set`. This has different
   * variants on which solutions to swap based on the FieldSolveType.
   */
  void
  update(FieldSolveType field_solve_type,
         Types::Index   solve_block,
         Types::Index   variable_index = 0);

  /**
   * @brief Prepare for solution transfer
   */
  void
  prepare_for_solution_transfer()
  {
    solution_transfer.prepare_for_coarsening_and_refinement(*solutions);
  }

  /**
   * @brief Transfer solutions
   */
  void
  execute_solution_transfer()
  {
    solution_transfer.interpolate(*solutions);
  }

  /**
   * @brief Reinit the solution transfer objections.
   */
  void
  reinit_solution_transfer(MatrixFreeContainer<dim, number> &matrix_free_container);

private:
  SolveGroup solve_group;
  // TODO (fractalsbyx): do this better
  int oldest_saved = 0;

  /**
   * @brief Mapping from block index to global field index.
   */
  std::vector<unsigned int> block_to_global_index;

  /**
   * @brief Mapping from global field index to block index.
   */
  std::vector<unsigned int> global_to_block_index;

  std::vector<FieldInfo::TensorRank> field_ranks;

  /**
   * @brief Whether multigrid has been enabled.
   */
  bool has_multigrid = false;

  /**
   * @brief Global minimum level for multigrid.
   */
  unsigned int global_min_level;

  /**
   * @brief Multigrid information.
   */
  const MGInfo<dim> *mg_info;

  /**
   * @brief The solution vectors.
   */
  struct SolutionLevel
  {
    BlockVector                                                      solutions;
    BlockVector                                                      new_solutions;
    std::array<BlockVector, Numbers::max_saved_increments>           old_solutions;
    dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> matrix_free;
  };

  // TODO (fractalsbyx): Consider switching to dealii::MGLevelObject
  // Right now, I prefer using the relative level.
  std::vector<SolutionLevel> solution_levels;

  /**
   * @brief Utility for solution transfer to different mesh (for AMR).
   */
  SolutionTransfer solution_transfer;
};

PRISMS_PF_END_NAMESPACE
