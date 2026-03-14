// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/dof_manager.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/system_wide.h>
#include <prismspf/core/timer.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
GroupSolutionHandler<dim, number>::GroupSolutionHandler(
  SolveGroup                          _solve_group,
  const std::vector<FieldAttributes> &_attributes_list)
  : solve_group(std::move(_solve_group))
{
  block_to_global_index.assign(solve_group.field_indices.begin(),
                               solve_group.field_indices.end());
  global_to_block_index =
    std::vector<Types::Index>(_attributes_list.size(), Numbers::invalid_index);
  for (unsigned int i = 0; i < block_to_global_index.size(); ++i)
    {
      global_to_block_index[block_to_global_index[i]] = i;
    }
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_solution_full_vector(unsigned int relative_level)
  -> BlockVector<number> &
{
  return solution_levels[relative_level].solutions;
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_solution_full_vector(
  unsigned int relative_level) const -> const BlockVector<number> &
{
  return solution_levels[relative_level].solutions;
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_solution_vector(unsigned int global_index,
                                                       unsigned int relative_level)
  -> SolutionVector<number> &
{
  return solution_levels[relative_level].solutions.block(
    global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_solution_vector(unsigned int global_index,
                                                       unsigned int relative_level) const
  -> const SolutionVector<number> &
{
  return solution_levels[relative_level].solutions.block(
    global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_old_solution_full_vector(
  unsigned int age,
  unsigned int relative_level) -> BlockVector<number> &
{
  return solution_levels[relative_level].old_solutions[age];
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_old_solution_full_vector(
  unsigned int age,
  unsigned int relative_level) const -> const BlockVector<number> &
{
  return solution_levels[relative_level].old_solutions[age];
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_old_solution_vector(unsigned int age,
                                                           unsigned int global_index,
                                                           unsigned int relative_level)
  -> SolutionVector<number> &
{
  return solution_levels[relative_level].old_solutions[age].block(
    global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_old_solution_vector(
  unsigned int age,
  unsigned int global_index,
  unsigned int relative_level) const -> const SolutionVector<number> &
{
  return solution_levels[relative_level].old_solutions[age].block(
    global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
SolutionLevel<dim, number> &
GroupSolutionHandler<dim, number>::get_solution_level(unsigned int relative_level)
{
  return solution_levels[relative_level];
}

template <unsigned int dim, typename number>
const SolutionLevel<dim, number> &
GroupSolutionHandler<dim, number>::get_solution_level(unsigned int relative_level) const
{
  return solution_levels[relative_level];
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_matrix_free(unsigned int relative_level)
  -> MatrixFree<dim, number> &
{
  return solution_levels[relative_level].matrix_free;
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_matrix_free(unsigned int relative_level) const
  -> const MatrixFree<dim, number> &
{
  return solution_levels[relative_level].matrix_free;
}

template <unsigned int dim, typename number>
unsigned int
GroupSolutionHandler<dim, number>::get_block_index(unsigned int global_index) const
{
  return global_to_block_index[global_index];
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_solve_group() const -> const SolveGroup &
{
  return solve_group;
}

template <unsigned int dim, typename number>
const std::vector<unsigned int> &
GroupSolutionHandler<dim, number>::get_global_to_block_index() const
{
  return global_to_block_index;
}

// TODO (fractalsbyx): This might all need to go in reinit(). Check if dof_handler and
// constraint ptrs change.
template <unsigned int dim, typename number>
template <unsigned int degree>
void
GroupSolutionHandler<dim, number>::init(
  const DoFManager<dim, degree>                &dof_manager,
  const ConstraintManager<dim, degree, number> &constraint_manager,
  unsigned int                                  num_old_saved)
{
  ConditionalOStreams::pout_base()
    << "Initializing solution set for solver " << solve_group.id << "...\n"
    << std::flush;
  Timer::start_section("initialize solution set");

  // TODO: figure out more consistent way of passing num_levels
  const unsigned int num_levels = dof_manager.get_dof_handlers().size();
  solution_levels.resize(num_levels);

  // Initialize matrixfree objects
  for (unsigned int relative_level = 0; relative_level < solution_levels.size();
       ++relative_level)
    {
      SolutionLevel<dim, number> &solution_level = solution_levels[relative_level];
      solution_level.old_solutions.resize(num_old_saved);
    }
  // Initialize solution vectors
  reinit(dof_manager, constraint_manager);
  // Initialize solution transfer
  init_solution_transfer();

  Timer::end_section("initialize solution set");
}

template <unsigned int dim, typename number>
template <unsigned int degree>
void
GroupSolutionHandler<dim, number>::reinit(
  const DoFManager<dim, degree>                &dof_manager,
  const ConstraintManager<dim, degree, number> &constraint_manager)
{
  // Initialize matrixfree objects
  for (unsigned int relative_level = 0; relative_level < solution_levels.size();
       ++relative_level)
    {
      solution_levels[relative_level].matrix_free.reinit(
        SystemWide<dim, degree>::mapping,
        dof_manager.get_field_dof_handlers(solve_group.field_indices, relative_level),
        constraint_manager.get_constraints(solve_group.field_indices, relative_level),
        dealii::QGaussLobatto<1>(degree + 1) // should dim really be 1?
      );
    }
  for (auto &solution_level : solution_levels)
    {
      BlockVector<number>              &solutions     = solution_level.solutions;
      std::vector<BlockVector<number>> &old_solutions = solution_level.old_solutions;
      MatrixFree<dim, number>          &matrix_free   = solution_level.matrix_free;

      // These partitioners basically just provide the number of elements in a distributed
      // way
      std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>>
        partitioners;
      partitioners.reserve(solve_group.field_indices.size());
      for (unsigned int block_index = 0; block_index < solve_group.field_indices.size();
           ++block_index)
        {
          partitioners.push_back(matrix_free.get_vector_partitioner(block_index));
        }
      // TODO (fractalsbyx): Check that the default MPI communicator is correct here
      solutions.reinit(partitioners);
      solutions.collect_sizes();
      for (BlockVector<number> &old_solution : old_solutions)
        {
          old_solution.reinit(partitioners);
          old_solution.collect_sizes();
        }
    }
}

// TODO (fractalsbyx): Check if this is necessary to repeat. Check if dof_handler ptrs
// change.
template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::init_solution_transfer()
{
  unsigned int num_blocks = solve_group.field_indices.size();
  block_solution_transfer.clear();
  block_solution_transfer.reserve(num_blocks);
  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
    {
      block_solution_transfer.emplace_back(
        get_matrix_free().get_dof_handler(block_index));
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::prepare_for_solution_transfer()
{
  unsigned int                      num_blocks    = solve_group.field_indices.size();
  const SolutionLevel<dim, number> &top_solutions = solution_levels[0];
  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
    {
      std::vector<const SolutionVector<number> *> fields_at_ages;
      fields_at_ages.reserve(1 + top_solutions.old_solutions.size());
      fields_at_ages.push_back(&(top_solutions.solutions.block(block_index)));
      for (const BlockVector<number> &old_solution : top_solutions.old_solutions)
        {
          fields_at_ages.push_back(&(old_solution.block(block_index)));
        }
      block_solution_transfer[block_index].prepare_for_coarsening_and_refinement(
        fields_at_ages);
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::execute_solution_transfer()
{
  // implementation will have identical structure to `prepare_for_solution_transfer()`
  unsigned int                num_blocks    = solve_group.field_indices.size();
  SolutionLevel<dim, number> &top_solutions = solution_levels[0];
  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
    {
      std::vector<SolutionVector<number> *> fields_at_ages;
      fields_at_ages.reserve(1 + top_solutions.old_solutions.size());
      fields_at_ages.push_back(&(top_solutions.solutions.block(block_index)));
      for (BlockVector<number> &old_solution : top_solutions.old_solutions)
        {
          fields_at_ages.push_back(&(old_solution.block(block_index)));
        }
      block_solution_transfer[block_index].interpolate(fields_at_ages);
    }
  update_ghosts(0);
  apply_constraints(0);
}

// TODO (fractalsbyx): Check if this is necessary for all solutions
template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::update_ghosts(unsigned int relative_level) const
{
  solution_levels[relative_level].solutions.update_ghost_values();
  for (const BlockVector<number> &old_solution :
       solution_levels[relative_level].old_solutions)
    {
      old_solution.update_ghost_values();
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::zero_out_ghosts(unsigned int relative_level) const
{
  solution_levels[relative_level].solutions.zero_out_ghost_values();
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::apply_constraints(unsigned int relative_level)
{
  BlockVector<number>     &solutions   = solution_levels[relative_level].solutions;
  MatrixFree<dim, number> &matrix_free = solution_levels[relative_level].matrix_free;

  unsigned int num_blocks = solve_group.field_indices.size();
  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
    {
      matrix_free.get_affine_constraints(block_index)
        .distribute(solutions.block(block_index));
    }
}

// TODO (fractalsbyx): Replace into initial_conditions module
template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::apply_initial_condition_for_old_fields()
{
  for (auto &solution_level : solution_levels)
    {
      BlockVector<number>              &solutions     = solution_level.solutions;
      std::vector<BlockVector<number>> &old_solutions = solution_level.old_solutions;

      for (BlockVector<number> &old_solution : old_solutions)
        {
          old_solution = solutions;
        }
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::update(unsigned int relative_level)
{
  // 1. TODO: propogate solutions to coarser levels (relative level needn't be an arg)

  // 3. bubble-swap method. bubble the discarded solution up to 'solution'
  SolutionLevel<dim, number> &solution_level = solution_levels[relative_level];
  for (int age = solution_level.old_solutions.size() - 1; age >= 0; --age)
    {
      if (age > 0)
        {
          solution_level.old_solutions[age].swap(solution_level.old_solutions[age - 1]);
        }
      else
        {
          solution_level.old_solutions[age].swap(solution_level.solutions);
        }
    }
}

#include "core/group_solution_handler.inst"

PRISMS_PF_END_NAMESPACE
