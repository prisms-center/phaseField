// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/memory_space.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <map>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Struct to hold the solution vectors.
 *
 * In PRISMS-PF we have several different allowable solves. At the core, each solve
 * updates some number of solution vectors using some number of old solutions.
 *
 * For explicit updates, we can represent it using a summation of old solutions with
 * arbitrary coefficients.
 *
 * \f[
 *  u^{n} = \sum_{j=1}^m a_j u^{n-j} = a_1 u^{n-1} + a_2 u^{n-2} ... a_j u^{n-j}
 * \f]
 *
 * \note \f$ m \f$ is finite and determined at compile time through
 * `prismspf::Numbers::max_saved_solutions`.
 *
 * For implicit updates, this becomes a little more complicated due to the use of
 * Newton's method, where \f$ u^n = u^{n-1} + \alpha^n \Delta u^n \f$. Here \f$ \alpha^n
 * \f$ is a damping parameter and \f$ \Delta u^n \f$ is the solution of the Newton step.
 * For each Newton step we solve for \f$ \Delta u^n \f$ and add it to the estimate of \f$
 * u^n \f$, iterating until we reach a satisfactory convergence.
 *
 * \f[
 *  u^{n} = a_0 \alpha^n \Delta u^n + \sum_{j=1}^m a_j u^{n-j} =  a_0 \alpha^n \Delta u^n
 * + a_1 u^{n-1} + a_2 u^{n-2} ... a_j u^{n-j}
 * \f]
 */
template <unsigned int dim, typename number>
struct SolutionBlock
{
  using BlockVector    = dealii::LinearAlgebra::distributed::BlockVector<number>;
  using SolutionVector = BlockVector::BlockType;

  /**
   * @brief Solutions being solved for.
   */
  BlockVector solutions;

  /**
   * @brief Old solution states.
   */
  std::array<BlockVector, Numbers::max_saved_solutions> old_solutions;

  /**
   * @brief Newton step solutions for implicit solves.
   */
  BlockVector change_solutions;
};

/**
 * @brief Class that manages solution initialization and swapping with old solutions.
 */
template <unsigned int dim, typename number>
class SolutionHandler

{
public:
  using VectorType   = dealii::LinearAlgebra::distributed::Vector<number>;
  using MGVectorType = dealii::LinearAlgebra::distributed::Vector<float>;

  /**
   * @brief Constructor.
   */
  SolutionHandler(const std::map<unsigned int, VariableAttributes> &_attributes_list,
                  const MGInfo<dim>                                &_mg_info);

  /**
   * @brief Get the solution vector set. This contains all the normal fields and is
   * typically used for output.
   *
   * TODO (landinjm): Make const ptr?
   */
  [[nodiscard]] std::map<unsigned int, VectorType *>
  get_solution_vector() const;

  /**
   * @brief Get a solution vector of a given field index and dependency type.
   *
   * TODO (landinjm): Make const ptr?
   */
  [[nodiscard]] VectorType *
  get_solution_vector(unsigned int index, DependencyType dependency_type) const;

  /**
   * @brief Get the "new" solution vector set.
   *
   * TODO (landinjm): Make const ptr?
   */
  [[nodiscard]] std::map<unsigned int, VectorType *>
  get_new_solution_vector() const;

  /**
   * @brief Get the "new" solution vector of a given field index.
   *
   * TODO (landinjm): Make const ptr?
   */
  [[nodiscard]] VectorType *
  get_new_solution_vector(unsigned int index) const;

  /**
   * @brief Get the mg solution vector set at a given level.
   */
  [[nodiscard]] std::vector<MGVectorType *>
  get_mg_solution_vector(unsigned int level) const;

  /**
   * @brief Get the mg solution vector set at a given level and index;
   */
  [[nodiscard]] MGVectorType *
  get_mg_solution_vector(unsigned int level, unsigned int index) const;

  /**
   * @brief Initialize the solution set.
   */
  void
  init(MatrixFreeContainer<dim, number> &matrix_free_container);

  /**
   * @brief Reinitialize the solution set.
   */
  void
  reinit(MatrixFreeContainer<dim, number> &matrix_free_container);

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
  apply_constraints(unsigned int                             index,
                    const dealii::AffineConstraints<number> &constraints);

  /**
   * @brief Apply initial condition to the old fields. For now, this simply copies the
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
    Assert(solution_transfer_set.size() == solution_set.size(),
           dealii::ExcInternalError());
    for (auto &[pair, solution] : solution_set)
      {
        auto &transfer = solution_transfer_set.at(pair);
        Assert(transfer, dealii::ExcInternalError());
        transfer->prepare_for_coarsening_and_refinement(*solution);
      }
  }

  /**
   * @brief Transfer solutions
   */
  void
  execute_solution_transfer()
  {
    Assert(solution_transfer_set.size() == solution_set.size(),
           dealii::ExcInternalError());
    for (auto &[pair, solution] : solution_set)
      {
        auto &transfer = solution_transfer_set.at(pair);
        Assert(transfer, dealii::ExcInternalError());
        transfer->interpolate(*solution);
      }
  }

  /**
   * @brief Free solution transfer objects.
   */
  void
  free_solution_transfer()
  {
    for (auto &[pair, ptr] : solution_transfer_set)
      {
        ptr.reset();
      }
    solution_transfer_set.clear();
  }

  /**
   * @brief Reinit the solution transfer objections.
   */
  void
  reinit_solution_transfer(MatrixFreeContainer<dim, number> &matrix_free_container);

private:
  /**
   * @brief The attribute list of the relevant variables.
   */
  const std::map<unsigned int, VariableAttributes> *attributes_list;

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
   * @brief The collection of solution vector at the current timestep. This includes
   * current values and old values.
   */
  std::map<std::pair<unsigned int, DependencyType>, std::unique_ptr<VectorType>>
    solution_set;

  /**
   * @brief Typedef for the solution transfer object.
   */
#if DEAL_II_VERSION_MAJOR >= 9 && DEAL_II_VERSION_MINOR >= 7
  using SolutionTransfer = dealii::SolutionTransfer<dim, VectorType>;
#else
  using SolutionTransfer =
    dealii::parallel::distributed::SolutionTransfer<dim, VectorType>;
#endif

  /**
   * @brief The collection of solution transfer objects at the current timestep.
   */
  std::map<std::pair<unsigned int, DependencyType>, std::unique_ptr<SolutionTransfer>>
    solution_transfer_set;

  /**
   * @brief The collection of new solution vectors at the current timestep. This is the
   * dst vector that is filled in the cell_loop. Unlike before, this only include the
   * current values which get updated in the `solution_set`.
   */
  std::map<unsigned int, std::unique_ptr<VectorType>> new_solution_set;

  /**
   * @brief The collection of solution vectors at the current timestep for the multigrid
   * hierarchy.
   */
  std::vector<std::vector<std::unique_ptr<MGVectorType>>> mg_solution_set;
};

PRISMS_PF_END_NAMESPACE
