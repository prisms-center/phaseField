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
  SolveBlock                            _solve_block,
  const std::vector<FieldAttributes>   &_attributes_list,
  const MatrixFreeManager<dim, number> &_matrix_free_manager)
  : solve_block(std::move(_solve_block))
  , matrix_free_manager(&_matrix_free_manager)
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
GroupSolutionHandler<dim, number>::get_solution_full_vector() -> BlockVector<number> &
{
  return primary_solutions.solutions;
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_solution_full_vector() const
  -> const BlockVector<number> &
{
  return primary_solutions.solutions;
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_solution_vector(unsigned int global_index)
  -> SolutionVector<number> &
{
  return primary_solutions.solutions.block(global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_solution_vector(unsigned int global_index) const
  -> const SolutionVector<number> &
{
  return primary_solutions.solutions.block(global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_old_solution_full_vector(unsigned int age)
  -> BlockVector<number> &
{
  return primary_solutions.old_solutions[age];
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_old_solution_full_vector(unsigned int age) const
  -> const BlockVector<number> &
{
  return primary_solutions.old_solutions[age];
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_old_solution_vector(unsigned int age,
                                                           unsigned int global_index)
  -> SolutionVector<number> &
{
  return primary_solutions.old_solutions[age].block(global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
auto
GroupSolutionHandler<dim, number>::get_old_solution_vector(
  unsigned int age,
  unsigned int global_index) const -> const SolutionVector<number> &
{
  return primary_solutions.old_solutions[age].block(global_to_block_index[global_index]);
}

template <unsigned int dim, typename number>
SolutionLevel<dim, number> &
GroupSolutionHandler<dim, number>::get_primary_solutions()
{
  return primary_solutions;
}

template <unsigned int dim, typename number>
const SolutionLevel<dim, number> &
GroupSolutionHandler<dim, number>::get_primary_solutions() const
{
  return primary_solutions;
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
  Assert(matrix_free_manager != nullptr, dealii::ExcNotInitialized());
  return matrix_free_manager->get_shared_matrix_free_levels();
}

template <unsigned int dim, typename number>
unsigned int
GroupSolutionHandler<dim, number>::get_block_index(unsigned int global_index) const
{
  return global_to_block_index[global_index];
}

template <unsigned int dim, typename number>
const SolveBlock &
GroupSolutionHandler<dim, number>::get_solve_block() const
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
GroupSolutionHandler<dim, number>::init(const NewDependencyExtents &extents)
{
  ConditionalOStreams::pout_base()
    << "Initializing solution set for solver " << solve_block.id << "...\n"
    << std::flush;
  Timer::start_section("Initialize solution set");

  primary_solutions.old_solutions.resize(extents.max_age);

  // Create solution levels
  const unsigned int num_levels = extents.max_age_per_level.size();
  solution_levels.resize(num_levels);

  // Make the correct number of old solution vectors
  for (unsigned int relative_level = 0; relative_level < solution_levels.size();
       ++relative_level)
    {
      SolutionLevel<dim, number> &solution_level = solution_levels[relative_level];
      solution_level.old_solutions.resize(extents.max_age_per_level[relative_level]);
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
  {
    BlockVector<number>              &solutions     = primary_solutions.solutions;
    std::vector<BlockVector<number>> &old_solutions = primary_solutions.old_solutions;

    // These partitioners basically just provide the number of elements in a distributed
    // way
    std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>> partitioners =
      matrix_free_manager->get_block_partitioners(solve_block.field_indices);
    solutions.reinit(partitioners);
    solutions.collect_sizes();
    for (BlockVector<number> &old_solution : old_solutions)
      {
        old_solution.reinit(partitioners);
        old_solution.collect_sizes();
      }
  }
  for (unsigned int relative_level = 0; relative_level < solution_levels.size();
       ++relative_level)
    {
      auto                             &solution_level = solution_levels[relative_level];
      BlockVector<number>              &solutions      = solution_level.solutions;
      std::vector<BlockVector<number>> &old_solutions  = solution_level.old_solutions;

      std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>>
        partitioners =
          matrix_free_manager->get_mg_block_partitioners(solve_block.field_indices,
                                                         relative_level);

      solutions.reinit(partitioners);
      solutions.collect_sizes();
      for (BlockVector<number> &old_solution : old_solutions)
        {
          old_solution.reinit(partitioners);
          old_solution.collect_sizes();
        }
    }
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
        matrix_free_manager->get_shared_matrix_free().get_dof_handler(field_index));
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::prepare_for_solution_transfer()
{
  unsigned int num_blocks = solve_block.field_indices.size();

  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
    {
      std::vector<const SolutionVector<number> *> fields_at_ages;
      fields_at_ages.reserve(1 + primary_solutions.old_solutions.size());
      fields_at_ages.push_back(&(primary_solutions.solutions.block(block_index)));
      for (const BlockVector<number> &old_solution : primary_solutions.old_solutions)
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
  unsigned int num_blocks = solve_block.field_indices.size();

  for (unsigned int block_index = 0; block_index < num_blocks; block_index++)
    {
      std::vector<SolutionVector<number> *> fields_at_ages;
      fields_at_ages.reserve(1 + primary_solutions.old_solutions.size());
      fields_at_ages.push_back(&(primary_solutions.solutions.block(block_index)));
      for (BlockVector<number> &old_solution : primary_solutions.old_solutions)
        {
          fields_at_ages.push_back(&(old_solution.block(block_index)));
        }
      block_solution_transfer[block_index].interpolate(fields_at_ages);
    }
  update_ghosts();
  apply_constraints_to_all();
}

// TODO (fractalsbyx): Check if this is necessary for all solutions
template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::update_ghosts() const
{
  primary_solutions.solutions.update_ghost_values();
  for (const BlockVector<number> &old_solution : primary_solutions.old_solutions)
    {
      old_solution.update_ghost_values();
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::zero_out_ghosts() const
{
  primary_solutions.solutions.zero_out_ghost_values();
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::apply_constraints_to_all()
{
  apply_constraints();
  std::vector<BlockVector<number>> &old_solutions = primary_solutions.old_solutions;
  const MatrixFree<dim, number>    &matrix_free =
    matrix_free_manager->get_shared_matrix_free();
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
GroupSolutionHandler<dim, number>::apply_constraints()
{
  apply_constraints(primary_solutions.solutions);
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::apply_constraints(BlockVector<number> &solution_vector)
{
  const MatrixFree<dim, number> &matrix_free =
    matrix_free_manager->get_shared_matrix_free();
  unsigned int num_blocks = solve_block.field_indices.size();
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
  BlockVector<number>              &solutions     = primary_solutions.solutions;
  std::vector<BlockVector<number>> &old_solutions = primary_solutions.old_solutions;
  for (BlockVector<number> &old_solution : old_solutions)
    {
      old_solution = solutions;
    }
}

// todo: refactor into solutionlevel
template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::update()
{
  // TODO: propagate solutions to coarser levels (relative level needn't be an arg)

  // bubble-swap method. bubble the discarded solution up to 'solution'
  for (int age = primary_solutions.old_solutions.size() - 1; age >= 0; --age)
    {
      if (age > 0)
        {
          primary_solutions.old_solutions[age].swap(
            primary_solutions.old_solutions[age - 1]);
        }
      else
        {
          primary_solutions.old_solutions[age].swap(primary_solutions.solutions);
        }
    }
  for (auto &solution_level : solution_levels)
    {
      for (int age = solution_level.old_solutions.size() - 1; age >= 0; --age)
        {
          if (age > 0)
            {
              solution_level.old_solutions[age].swap(
                solution_level.old_solutions[age - 1]);
            }
          else
            {
              solution_level.old_solutions[age].swap(solution_level.solutions);
            }
        }
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::print_solution_full_vector(std::ostream &out) const
{
  // Print the solutions with 12 digits of precision
  get_solution_full_vector().print(out, 12);
}

template <unsigned int dim, typename number>
unsigned int
GroupSolutionHandler<dim, number>::num_levels()
{
  return solution_levels.size();
}

#include "core/group_solution_handler.inst"

PRISMS_PF_END_NAMESPACE
