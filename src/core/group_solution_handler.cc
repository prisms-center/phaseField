// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

// #include <deal.II/base/exceptions.h>
// #include <deal.II/base/mg_level_object.h>
// #include <deal.II/distributed/solution_transfer.h>
// #include <deal.II/lac/affine_constraints.h>
// #include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/group_solution_handler.h>
// #include <prismspf/core/matrix_free_handler.h>
// #include <prismspf/core/multigrid_info.h>
// #include <prismspf/core/solution_handler.h> //
#include <prismspf/core/timer.h>
// #include <prismspf/core/type_enums.h>
// #include <prismspf/core/types.h>
// #include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
GroupSolutionHandler<dim, number>::GroupSolutionHandler(
  const SolveGroup                   &_solve_group,
  const std::vector<FieldAttributes> &_attributes_list,
  const MGInfo<dim>                  &_mg_info)
  : mg_info(&_mg_info)
{
  block_to_global_index.assign(_solve_group.field_indices.begin(),
                               _solve_group.field_indices.end());
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
void GroupSolutionHandler<dim, number>::init(/* args */)
{
  ConditionalOStreams::pout_base()
    << "Initializing solution set for solver " << solve_group.id << "...\n"
    << std::flush;
  Timer::start_section("reinitialize solution set");

  // Initialize matrixfree objects
  // TODO
  for (auto &solution_level : solution_levels)
    {
      dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &matrix_free =
        solution_level.matrix_free;
      matrix_free.reinit(/* args */);
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
      dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &matrix_free =
        solution_level.matrix_free;

      // These partitioners basically just provide the number of elements in a distributed
      // way
      std::vector<std::shared_ptr<const dealii::Utilities::MPI::Partitioner>>
        partitioners;
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

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::reinit_solution_transfer(
  MatrixFreeContainer<dim, number> &matrix_free_container)
{
  // Create all entries
  for (const auto &[index, variable] : *attributes_list)
    {
      // Add the variable if it doesn't already exist
      if (!solution_transfer_set.contains(std::make_pair(index, DependencyType::Normal)))
        {
          solution_transfer_set[std::make_pair(index, DependencyType::Normal)] =
            std::make_unique<SolutionTransfer>(
              matrix_free_container.get_matrix_free()->get_dof_handler(index));
        }

      // Add dependencies if they don't exist
      Types::Index field_index = 0;
      for (const auto &dependency_set : variable.get_eval_flag_set_rhs())
        {
          Types::Index dep_index = 0;
          for (const auto &value : dependency_set)
            {
              if (value == dealii::EvaluationFlags::EvaluationFlags::nothing)
                {
                  dep_index++;
                  continue;
                }

              solution_transfer_set.try_emplace(
                std::make_pair(field_index, static_cast<DependencyType>(dep_index)),
                std::make_unique<SolutionTransfer>(
                  matrix_free_container.get_matrix_free()->get_dof_handler(field_index)));

              dep_index++;
            }

          field_index++;
        }
      field_index = 0;
      for (const auto &dependency_set : variable.get_eval_flag_set_lhs())
        {
          Types::Index dep_index = 0;
          for (const auto &value : dependency_set)
            {
              if (value == dealii::EvaluationFlags::EvaluationFlags::nothing)
                {
                  dep_index++;
                  continue;
                }

              solution_transfer_set.try_emplace(
                std::make_pair(field_index, static_cast<DependencyType>(dep_index)),
                std::make_unique<SolutionTransfer>(
                  matrix_free_container.get_matrix_free()->get_dof_handler(field_index)));

              dep_index++;
            }

          field_index++;
        }
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::update_ghosts() const
{
  for (const auto &[pair, solution] : solution_set)
    {
      solution->update_ghost_values();
    }
  for (const auto &index_vector : mg_solution_set)
    {
      for (const auto &solution : index_vector)
        {
          solution->update_ghost_values();
        }
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::zero_out_ghosts() const
{
  for (const auto &[index, solution] : new_solution_set)
    {
      solution->zero_out_ghost_values();
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::apply_constraints(
  unsigned int                             index,
  const dealii::AffineConstraints<number> &constraints)
{
  for (auto &[pair, vector] : solution_set)
    {
      if (pair.first != index)
        {
          continue;
        }
      constraints.distribute(*vector);
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::apply_initial_condition_for_old_fields()
{
  for (auto &[pair, vector] : solution_set)
    {
      if (pair.second == DependencyType::Normal)
        {
          continue;
        }
      *(get_solution_vector(pair.first, pair.second)) =
        *(get_solution_vector(pair.first, DependencyType::Normal));
    }
}

template <unsigned int dim, typename number>
void
GroupSolutionHandler<dim, number>::update(FieldSolveType field_solve_type,
                                          Types::Index   solve_block,
                                          Types::Index   variable_index)
{
  // Helper function to swap vectors for all dependency types
  auto swap_all_dependency_vectors = [this](Types::Index index, auto &new_vector)
  {
    // Always swap the Normal dependency
    new_vector->swap(*(solution_set.at(std::make_pair(index, DependencyType::Normal))));

    // Swap old dependency types if they exist
    const std::array<DependencyType, 4> old_types = {
      {DependencyType::OldOne,
       DependencyType::OldTwo,
       DependencyType::OldThree,
       DependencyType::OldFour}
    };

    for (const auto &dep_type : old_types)
      {
        if (solution_set.contains(std::make_pair(index, dep_type)))
          {
            new_vector->swap(*(solution_set.at(std::make_pair(index, dep_type))));
          }
      }
  };

  // Loop through the solutions and swap them
  for (auto &[index, new_vector] : new_solution_set)
    {
      const auto &attr_field_type = attributes_list->at(index).get_field_solve_type();

      // Skip if the solve block is wrong
      if (attributes_list->at(index).get_solve_block() != solve_block)
        {
          continue;
        }

      // Skip if field solve types don't match
      if (attr_field_type != field_solve_type)
        {
          continue;
        }

      switch (field_solve_type)
        {
          case FieldSolveType::ExplicitConstant:
            // For ExplicitConstant we don't do anything
            break;
          case FieldSolveType::Explicit:
          case FieldSolveType::ExplicitPostprocess:
            // For ExplicitPostprocess we only swap Normal, but the helper function will
            // do that first and ignore the rest since they should exist
            swap_all_dependency_vectors(index, new_vector);
            break;
          case FieldSolveType::NonexplicitLinear:
          case FieldSolveType::NonexplicitAuxiliary:
          case FieldSolveType::NonexplicitSelfnonlinear:
          case FieldSolveType::NonexplicitCononlinear:
            // For Nonexplicit types we swap all dependency types if the index matches
            if (variable_index == index)
              {
                swap_all_dependency_vectors(index, new_vector);
              }
            break;
          default:
            AssertThrow(false, dealii::ExcMessage("Invalid FieldSolveType"));
        }
    }
}

#include "core/solution_handler.inst"

PRISMS_PF_END_NAMESPACE
