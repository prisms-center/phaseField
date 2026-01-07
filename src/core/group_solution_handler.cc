// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/timer.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
GroupSolutionHandler<dim, number>::GroupSolutionHandler(
  const SolveGroup                   &_solve_group,
  const std::vector<FieldAttributes> &_attributes_list)
{
  block_to_global_index.assign(_solve_group.field_indices.begin(),
                               _solve_group.field_indices.end());
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
  -> BlockVector &
{
  return solution_levels[relative_level].solutions;
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_solution_vector(unsigned int global_index,
                                                       unsigned int relative_level)
  -> SolutionVector &
{
  return solution_levels[relative_level].solutions.block(
    global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_old_solution_full_vector(
  unsigned int age,
  unsigned int relative_level) -> BlockVector &
{
  return solution_levels[relative_level].old_solutions[age];
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_old_solution_vector(unsigned int age,
                                                           unsigned int global_index,
                                                           unsigned int relative_level)
  -> SolutionVector &
{
  return solution_levels[relative_level].old_solutions[age].block(
    global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_new_solution_full_vector(
  unsigned int relative_level) -> BlockVector &
{
  return solution_levels[relative_level].new_solutions;
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_new_solution_vector(unsigned int global_index,
                                                           unsigned int relative_level)
  -> SolutionVector &
{
  return solution_levels[relative_level].new_solutions.block(
    global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
SolutionLevel<dim, number> &
GroupSolutionHandler<dim, number>::get_solution_level(unsigned int relative_level)
{
  return solution_levels[relative_level];
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_matrix_free(unsigned int relative_level)
  -> MatrixFree &
{
  return solution_levels[relative_level].matrix_free;
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_solution_full_vector(
  unsigned int relative_level) const -> const BlockVector &
{
  return solution_levels[relative_level].solutions;
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_solution_vector(unsigned int global_index,
                                                       unsigned int relative_level) const
  -> const SolutionVector &
{
  return solution_levels[relative_level].solutions.block(
    global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_old_solution_full_vector(
  unsigned int age,
  unsigned int relative_level) const -> const BlockVector &
{
  return solution_levels[relative_level].old_solutions[age];
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_old_solution_vector(
  unsigned int age,
  unsigned int global_index,
  unsigned int relative_level) const -> const SolutionVector &
{
  return solution_levels[relative_level].old_solutions[age].block(
    global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_new_solution_full_vector(
  unsigned int relative_level) const -> const BlockVector &
{
  return solution_levels[relative_level].new_solutions;
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_new_solution_vector(
  unsigned int global_index,
  unsigned int relative_level) const -> const SolutionVector &
{
  return solution_levels[relative_level].new_solutions.block(
    global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
const SolutionLevel<dim, number> &
GroupSolutionHandler<dim, number>::get_solution_level(unsigned int relative_level) const
{
  return solution_levels[relative_level];
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_matrix_free(unsigned int relative_level) const
  -> const MatrixFree &
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
  const dealii::Mapping<dim>                   &mapping,
  const DofManager<dim>                        &dof_manager,
  const ConstraintManager<dim, degree, number> &constraint_manager,
  const dealii::Quadrature<dim>                &quad)
{
  ConditionalOStreams::pout_base()
    << "Initializing solution set for solver " << solve_group.id << "...\n"
    << std::flush;
  Timer::start_section("reinitialize solution set");

  // Initialize matrixfree objects
  // TODO
  for (unsigned int relative_level = 0; relative_level < solution_levels.size();
       ++relative_level)
    {
      SolutionLevel<dim, number> &solution_level = solution_levels[relative_level];
      MatrixFree                 &matrix_free    = solution_level.matrix_free;
      matrix_free.reinit(
        mapping,
        dof_manager.get_dof_handlers(solve_group.field_indices, relative_level),
        constraint_manager.get_constraints(solve_group.field_indices, relative_level),
        quad);
    }
  // Initialize solution vectors
  reinit();

  Timer::end_section("reinitialize solution set");
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::reinit()
{
  for (auto &solution_level : solution_levels)
    {
      BlockVector &solutions     = solution_level.solutions;
      BlockVector &new_solutions = solution_level.new_solutions;
      std::array<BlockVector, Numbers::max_saved_increments> &old_solutions =
        solution_level.old_solutions;
      MatrixFree &matrix_free = solution_level.matrix_free;

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
      new_solutions.reinit(partitioners);
      for (unsigned int i = 0; i < Numbers::max_saved_increments; ++i)
        {
          if (i < oldest_saved)
            {
              old_solutions[i].reinit(partitioners);
            }
        }
      // TODO (fractalsbyx): Check if ghosts need to be updated here
    }
}

// TODO (fractalsbyx): Check if this is necessary to repeat. Check if dof_handler ptrs
// change.
template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::init_solution_transfer()
{
  const dealii::DoFHandler<dim> &dof_handler =
    solution_levels[0].matrix_free.get_dof_handler();
  solution_transfer = SolutionTransfer(dof_handler);
  for (unsigned int i = 0; i < oldest_saved; ++i)
    {
      old_solution_transfer[i] = SolutionTransfer(dof_handler);
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::prepare_for_solution_transfer()
{
  solution_transfer.prepare_for_coarsening_and_refinement(solution_levels[0].solutions);
  for (unsigned int i = 0; i < oldest_saved; ++i)
    {
      old_solution_transfer[i].prepare_for_coarsening_and_refinement(
        solution_levels[0].old_solutions[i]);
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::execute_solution_transfer()
{
  solution_transfer.interpolate(solution_levels[0].solutions);
  for (unsigned int i = 0; i < oldest_saved; ++i)
    {
      old_solution_transfer[i].interpolate(solution_levels[0].old_solutions[i]);
    }
}

// TODO (fractalsbyx): Check if this is necessary for all solutions
template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::update_ghosts(unsigned int relative_level) const
{
  solution_levels[relative_level].solutions.update_ghost_values();
  solution_levels[relative_level].new_solutions.update_ghost_values();
  for (unsigned int i = 0; i < oldest_saved; ++i)
    {
      solution_levels[relative_level].old_solutions[i].update_ghost_values();
    }
}

// TODO (fractalsbyx): Check if this is necessary for all solutions
template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::zero_out_ghosts(unsigned int relative_level) const
{
  solution_levels[relative_level].solutions.zero_out_ghost_values();
  solution_levels[relative_level].new_solutions.zero_out_ghost_values();
  for (unsigned int i = 0; i < oldest_saved; ++i)
    {
      solution_levels[relative_level].old_solutions[i].zero_out_ghost_values();
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::apply_constraints(unsigned int relative_level)
{
  BlockVector &solutions     = solution_levels[relative_level].solutions;
  BlockVector &new_solutions = solution_levels[relative_level].new_solutions;
  MatrixFree  &matrix_free   = solution_levels[relative_level].matrix_free;

  for (unsigned int i = 0; i < matrix_free.get_affine_constraints().size(); ++i)
    {
      matrix_free.get_affine_constraints()[i].distribute(solutions.block(i));
      matrix_free.get_affine_constraints()[i].distribute(new_solutions.block(i));
      // TODO: Check if this needs to be done to both or just one of the above
    }
}

// TODO (fractalsbyx): Replace into initial_conditions module
template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::apply_initial_condition_for_old_fields()
{
  for (auto &solution_level : solution_levels)
    {
      BlockVector &solutions = solution_level.solutions;
      std::array<BlockVector, Numbers::max_saved_increments> &old_solutions =
        solution_level.old_solutions;

      for (unsigned int i = 0; i < oldest_saved; ++i)
        {
          old_solutions[i] = solutions;
        }
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::update(unsigned int relative_level)
{
  // bubble-swap method. bubble the discarded solution up to 'new_solution'
  SolutionLevel<dim, number> &solution_level = solution_levels[relative_level];
  for (int age = oldest_saved - 1; age >= 0; --age)
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
  solution_level.solutions.swap(solution_level.new_solutions);
}

// #include "core/group_solution_handler.inst"

PRISMS_PF_END_NAMESPACE
