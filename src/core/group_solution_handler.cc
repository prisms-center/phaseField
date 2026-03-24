// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/dof_manager.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/group_solution_handler.h>
#include <prismspf/core/solve_block.h>
#include <prismspf/core/system_wide.h>
#include <prismspf/core/timer.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
GroupSolutionHandler<dim, number>::GroupSolutionHandler(
  SolveBlock                                  _solve_block,
  const std::vector<FieldAttributes>         &_attributes_list,
  const std::vector<MatrixFree<dim, number>> &_matrix_free_levels)
  : solve_block(std::move(_solve_block))
  , matrix_free_levels(&_matrix_free_levels)
  , mg_transfer(dealii::MGTransferMF<dim, number>())
{
  block_to_global_index.assign(solve_block.field_indices.begin(),
                               solve_block.field_indices.end());
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
const std::vector<MatrixFree<dim, number>> &
GroupSolutionHandler<dim, number>::get_matrix_free_levels() const
{
  Assert(matrix_free_levels != nullptr, dealii::ExcNotInitialized());
  return *matrix_free_levels;
}

template <unsigned int dim, typename number>
unsigned int
GroupSolutionHandler<dim, number>::get_block_index(unsigned int global_index) const
{
  return global_to_block_index[global_index];
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_solve_block() const -> const SolveBlock &
{
  return solve_block;
}

template <unsigned int dim, typename number>
const std::vector<unsigned int> &
GroupSolutionHandler<dim, number>::get_global_to_block_index() const
{
  return global_to_block_index;
}

template <unsigned int dim, typename number>
const std::vector<unsigned int> &
GroupSolutionHandler<dim, number>::get_block_to_global_index() const
{
  return block_to_global_index;
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::init(
  const std::vector<unsigned int> &max_age_per_level)
{
  ConditionalOStreams::pout_base()
    << "Initializing solution set for solver " << solve_block.id << "...\n"
    << std::flush;
  Timer::start_section("Initialize solution set");

  // Create solution levels
  const unsigned int num_levels = max_age_per_level.size();
  solution_levels.resize(num_levels);

  // Make the correct number of old solution vectors
  for (unsigned int relative_level = 0; relative_level < solution_levels.size();
       ++relative_level)
    {
      SolutionLevel<dim, number> &solution_level = solution_levels[relative_level];
      solution_level.old_solutions.resize(max_age_per_level[relative_level]);
    }

  // Initialize solution vectors
  reinit();
  // Initialize solution transfer
  init_solution_transfer();

  Timer::end_section("Initialize solution set");
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::reinit()
{
  for (unsigned int relative_level = 0; relative_level < solution_levels.size();
       ++relative_level)
    {
      auto                             &solution_level = solution_levels[relative_level];
      BlockVector<number>              &solutions      = solution_level.solutions;
      std::vector<BlockVector<number>> &old_solutions  = solution_level.old_solutions;
      const MatrixFree<dim, number>    &matrix_free =
        get_matrix_free_levels()[relative_level];

      // These partitioners basically just provide the number of elements in a distributed
      // way
      std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>>
        partitioners;
      partitioners.reserve(solve_block.field_indices.size());
      for (unsigned int field_index : solve_block.field_indices)
        {
          partitioners.push_back(matrix_free.get_vector_partitioner(field_index));
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
  // mg_transfer.build(dof_manager.get_field_dof_handlers(solve_block.field_indices, 0));
  // todo
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::init_solution_transfer()
{
  unsigned int num_blocks = solve_block.field_indices.size();
  block_solution_transfer.clear();
  block_solution_transfer.reserve(num_blocks);
  for (unsigned int field_index : solve_block.field_indices)
    {
      block_solution_transfer.emplace_back(
        get_matrix_free_levels()[0].get_dof_handler(field_index));
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::prepare_for_solution_transfer()
{
  unsigned int                      num_blocks    = solve_block.field_indices.size();
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
  unsigned int                num_blocks    = solve_block.field_indices.size();
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
  apply_constraints_to_all(0);
}

template <unsigned int dim, typename number>
template <unsigned int degree>
void
GroupSolutionHandler<dim, number>::mg_transfer_down(
  const DoFManager<dim, degree> &dof_manager,
  unsigned int                   finest_level,
  bool                           transfer_old_solutions)
{
  dealii::MGLevelObject<BlockVector<number>> temp_mg_solutions(1 + finest_level -
                                                                 solution_levels.size(),
                                                               finest_level);

  const auto dof_handlers =
    dof_manager.get_block_dof_handlers(solve_block.field_indices, 0);

  // transfer regular solutions to mg levels
  mg_transfer.copy_to_mg(dof_handlers, temp_mg_solutions, solution_levels[0].solutions);
  // swap to actual mg solution vectors
  for (unsigned int relative_level = 0; relative_level < solution_levels.size();
       ++relative_level)
    {
      unsigned int level = finest_level - relative_level;
      solution_levels[relative_level].solutions.swap(temp_mg_solutions[level]);
    }
  if (!transfer_old_solutions)
    {
      return;
    }
  // transfer old solutions
  for (unsigned int age_index = 0; age_index < solution_levels[0].old_solutions.size();
       ++age_index)
    {
      mg_transfer.copy_to_mg(dof_handlers,
                             temp_mg_solutions,
                             solution_levels[0].old_solutions[age_index]);
      for (unsigned int relative_level = 0; relative_level < solution_levels.size();
           ++relative_level)
        {
          unsigned int level = finest_level - relative_level;
          if (age_index < solution_levels[relative_level].old_solutions.size())
            {
              solution_levels[relative_level].old_solutions[age_index].swap(
                temp_mg_solutions[level]);
            }
        }
    }
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
GroupSolutionHandler<dim, number>::apply_constraints_to_all(unsigned int relative_level)
{
  apply_constraints(0);
  std::vector<BlockVector<number>> &old_solutions =
    solution_levels[relative_level].old_solutions;
  const MatrixFree<dim, number> &matrix_free = get_matrix_free_levels()[relative_level];
  for (BlockVector<number> &solutions : old_solutions)
    {
      unsigned int num_blocks = solve_block.field_indices.size();
      for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
        {
          matrix_free.get_affine_constraints(block_to_global_index[block_index])
            .distribute(solutions.block(block_index));
        }
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::apply_constraints(unsigned int relative_level)
{
  BlockVector<number> &solutions = solution_levels[relative_level].solutions;
  apply_constraints(solutions, relative_level);
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::apply_constraints(BlockVector<number> &solution_vector,
                                                     unsigned int         relative_level)
{
  const MatrixFree<dim, number> &matrix_free = get_matrix_free_levels()[relative_level];
  unsigned int                   num_blocks  = solve_block.field_indices.size();
  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
    {
      matrix_free.get_affine_constraints(block_to_global_index[block_index])
        .distribute(solution_vector.block(block_index));
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
  // 1. TODO: propagate solutions to coarser levels (relative level needn't be an arg)

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

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::print_solution_full_vector(
  std::ostream &out,
  unsigned int  relative_level) const
{
  // Print the solutions with 12 digits of precision
  get_solution_full_vector(relative_level).print(out, 12);
}

#include "core/group_solution_handler.inst"

PRISMS_PF_END_NAMESPACE
