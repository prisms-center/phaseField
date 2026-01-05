// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/constraint_manager.h>
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief The solution vectors.
 */
template <unsigned int dim, typename number>
struct SolutionLevel
{
  using BlockVector    = dealii::LinearAlgebra::distributed::BlockVector<number>;
  using SolutionVector = BlockVector::BlockType;
  using MatrixFree     = dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>;
  BlockVector                                            solutions;
  BlockVector                                            new_solutions;
  std::array<BlockVector, Numbers::max_saved_increments> old_solutions;
  BlockVector                                            change_solutions;
  MatrixFree                                             matrix_free;
};

/**
 * @brief Class that manages solution initialization and swapping with old solutions.
 */
template <unsigned int dim, typename number>
class GroupSolutionHandler
{
public:
  using BlockVector    = SolutionLevel<dim, number>::BlockVector;
  using SolutionVector = SolutionLevel<dim, number>::SolutionVector;
  using MatrixFree     = SolutionLevel<dim, number>::MatrixFree;
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
                       const std::vector<FieldAttributes> &_attributes_list);

  // TODO (fractalsbyx): add 'const T& get() const' versions.
  /**
   * @brief Get the solution vector set. This contains all the normal fields and is
   * typically used for output.
   */
  [[nodiscard]] BlockVector &
  get_solution_full_vector(unsigned int relative_level = 0);

  /**
   * @brief Get the const solution vector set. This contains all the normal fields and is
   * typically used for output.
   */
  [[nodiscard]] const BlockVector &
  get_solution_full_vector(unsigned int relative_level = 0) const;

  /**
   * @brief Get a solution vector of a given field index.
   */
  [[nodiscard]] SolutionVector &
  get_solution_vector(unsigned int global_index, unsigned int relative_level = 0);

  /**
   * @brief Get a solution vector of a given field index.
   */
  [[nodiscard]] const SolutionVector &
  get_solution_vector(unsigned int global_index, unsigned int relative_level = 0) const;

  /**
   * @brief Get the old solution vector set at a given age.
   */
  [[nodiscard]] BlockVector &
  get_old_solution_full_vector(unsigned int age, unsigned int relative_level = 0);

  /**
   * @brief Get the old solution vector set at a given age.
   */
  [[nodiscard]] const BlockVector &
  get_old_solution_full_vector(unsigned int age, unsigned int relative_level = 0) const;

  /**
   * @brief Get a solution vector of a given field index at a given age.
   */
  [[nodiscard]] SolutionVector &
  get_old_solution_vector(unsigned int age,
                          unsigned int global_index,
                          unsigned int relative_level = 0);

  /**
   * @brief Get a solution vector of a given field index at a given age.
   */
  [[nodiscard]] const SolutionVector &
  get_old_solution_vector(unsigned int age,
                          unsigned int global_index,
                          unsigned int relative_level = 0) const;

  /**
   * @brief Get the "new" solution vector set.
   */
  [[nodiscard]] BlockVector &
  get_new_solution_full_vector(unsigned int relative_level = 0);

  /**
   * @brief Get the "new" solution vector set.
   */
  [[nodiscard]] const BlockVector &
  get_new_solution_full_vector(unsigned int relative_level = 0) const;

  /**
   * @brief Get the "new" solution vector of a given field index.
   */
  [[nodiscard]] SolutionVector &
  get_new_solution_vector(unsigned int index, unsigned int relative_level = 0);

  /**
   * @brief Get the "new" solution vector of a given field index.
   */
  [[nodiscard]] const SolutionVector &
  get_new_solution_vector(unsigned int index, unsigned int relative_level = 0) const;

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
   * @brief Get the matrix_free object at a level.
   */
  [[nodiscard]] MatrixFree &
  get_matrix_free(unsigned int relative_level = 0);

  /**
   * @brief Get the matrix_free object at a level.
   */
  [[nodiscard]] const MatrixFree &
  get_matrix_free(unsigned int relative_level = 0) const;

  /**
   * @brief Get the block index from the global index.
   */
  [[nodiscard]] unsigned int
  get_block_index(unsigned int global_index) const;

  /**
   * @brief Get the underlying solve group object.
   */
  [[nodiscard]] const SolveGroup &
  get_solve_group() const;

  /**
   * @brief Initialize the solution set.
   */
  template <unsigned int degree>
  void
  init(const dealii::Mapping<dim>                   &mapping,
       const DofManager<dim>                        &dof_manager,
       const ConstraintManager<dim, degree, number> &constraint_manager,
       const dealii::Quadrature<dim>                &quad);

  /**
   * @brief Reinitialize the solution set.
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
   * @brief Apply the given constraints to a solution vector of a given field index.
   */
  void
  apply_constraints(unsigned int relative_level = 0);

  /**
   * @brief Apply intial condition to the old fields. For now, this simply copies the
   * values in the normal field to the old.
   */
  void
  apply_initial_condition_for_old_fields();

  /**
   * @brief Update the `solutions` to the `new_solution` and propagate the old solutions.
   */
  void
  update(unsigned int relative_level = 0);

  /**
   * @brief Reinit the solution transfer objections.
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

private:
  /**
   * @brief Information about the solve group this handler is responsible for.
   */
  SolveGroup solve_group;

  // TODO (fractalsbyx): do this better
  /**
   * @brief Oldest saved age.
   */
  int oldest_saved = 0;

  /**
   * @brief Mapping from block index to global field index.
   */
  std::vector<unsigned int> block_to_global_index;

  /**
   * @brief Mapping from global field index to block index.
   */
  std::vector<unsigned int> global_to_block_index;

  // TODO (fractalsbyx): Consider switching to dealii::MGLevelObject
  // Right now, I prefer using the relative level.
  std::vector<SolutionLevel<dim, number>> solution_levels;

  /**
   * @brief Utility for solution transfer to different mesh (for AMR).
   */
  SolutionTransfer                                            solution_transfer;
  std::array<SolutionTransfer, Numbers::max_saved_increments> old_solution_transfer;
};

PRISMS_PF_END_NAMESPACE
