// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/concurrent_solver.h>
#include <prismspf/solvers/solver_base.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <functional>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
ConcurrentSolver<dim, degree, number>::ConcurrentSolver(
  const SolverContext<dim, degree, number> &_solver_context,
  const FieldSolveType                     &_field_solve_type,
  Types::Index                              _solve_priority)
  : SolverBase<dim, degree, number>(_solver_context, _field_solve_type, _solve_priority)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentSolver<dim, degree, number>::init()
{
  // Call the base class init
  this->SolverBase<dim, degree, number>::init();

  // If the FieldSolveType is constant or the solver is empty we can just return early.
  if (this->get_field_solve_type() == FieldSolveType::ExplicitConstant ||
      this->solver_is_empty())
    {
      return;
    }

  // Create the MatrixFreeOperator
  system_matrix =
    std::make_unique<typename SolverBase<dim, degree, number>::SystemMatrixType>(
      this->get_subset_attributes(),
      this->get_pde_operator(),
      this->get_solve_block());

  // Set up the user-implemented equations and create the residual vectors
  system_matrix->clear();
  system_matrix->initialize(this->get_matrix_free_container().get_matrix_free(),
                            this->get_element_volume_container().get_element_volume());

  // Grab some data from the VariableAttributes
  const Types::Index max_fields =
    this->get_subset_attributes().begin()->second.get_max_fields();
  const Types::Index max_dependency_types =
    this->get_subset_attributes().begin()->second.get_max_dependency_types();

  // Resize the global to local solution vector
  global_to_local_solution.resize(static_cast<unsigned long>(max_fields) *
                                    max_dependency_types,
                                  Numbers::invalid_index);

  // Create the subset of solution vectors and add the mapping to MatrixFreeOperator
  Types::Index dependency_index = 0;
  for (const auto &inner_dependency_set :
       this->get_subset_attributes().begin()->second.get_dependency_set_rhs())
    {
      Types::Index dependency_type = 0;
      for (const auto &field_type : inner_dependency_set)
        {
          // Skip if an invalid field type is found or the global_to_local_solution
          // already has an entry for this dependency index and dependency type
          if (field_type == FieldInfo::TensorRank::Undefined ||
              global_to_local_solution[(dependency_index * max_dependency_types) +
                                       dependency_type] != Numbers::invalid_index)
            {
              dependency_type++;
              continue;
            }

          solution_subset.push_back(this->get_solution_handler().get_solution_vector(
            dependency_index,
            static_cast<DependencyType>(dependency_type)));
          new_solution_subset.push_back(
            this->get_solution_handler().get_new_solution_vector(dependency_index));
          global_to_local_solution[(dependency_index * max_dependency_types) +
                                   dependency_type] =
            static_cast<unsigned int>(solution_subset.size()) - 1;

          dependency_type++;
        }

      dependency_index++;
    }
  system_matrix->add_global_to_local_mapping(global_to_local_solution);
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentSolver<dim, degree, number>::reinit()
{
  // Call the base class reinit
  this->SolverBase<dim, degree, number>::reinit();

  // If the FieldSolveType is constant or the solver is empty we can just return early.
  if (this->get_field_solve_type() == FieldSolveType::ExplicitConstant ||
      this->solver_is_empty())
    {
      return;
    }

  // Clear some objects
  global_to_local_solution.clear();
  solution_subset.clear();
  new_solution_subset.clear();

  // Set up the user-implemented equations and create the residual vectors
  system_matrix->clear();
  system_matrix->initialize(this->get_matrix_free_container().get_matrix_free(),
                            this->get_element_volume_container().get_element_volume());

  // Grab some data from the VariableAttributes
  const Types::Index max_fields =
    this->get_subset_attributes().begin()->second.get_max_fields();
  const Types::Index max_dependency_types =
    this->get_subset_attributes().begin()->second.get_max_dependency_types();

  // Resize the global to local solution vector
  global_to_local_solution.resize(static_cast<unsigned long>(max_fields) *
                                    max_dependency_types,
                                  Numbers::invalid_index);

  // Create the subset of solution vectors and add the mapping to MatrixFreeOperator
  Types::Index dependency_index = 0;
  for (const auto &inner_dependency_set :
       this->get_subset_attributes().begin()->second.get_dependency_set_rhs())
    {
      Types::Index dependency_type = 0;
      for (const auto &field_type : inner_dependency_set)
        {
          // Skip if an invalid field type is found or the global_to_local_solution
          // already has an entry for this dependency index and dependency type
          if (field_type == FieldInfo::TensorRank::Undefined ||
              global_to_local_solution[(dependency_index * max_dependency_types) +
                                       dependency_type] != Numbers::invalid_index)
            {
              dependency_type++;
              continue;
            }

          solution_subset.push_back(this->get_solution_handler().get_solution_vector(
            dependency_index,
            static_cast<DependencyType>(dependency_type)));
          new_solution_subset.push_back(
            this->get_solution_handler().get_new_solution_vector(dependency_index));
          global_to_local_solution[(dependency_index * max_dependency_types) +
                                   dependency_type] =
            static_cast<unsigned int>(solution_subset.size()) - 1;

          dependency_type++;
        }

      dependency_index++;
    }
  system_matrix->add_global_to_local_mapping(global_to_local_solution);
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentSolver<dim, degree, number>::solve()
{
  // Call the base class solve
  this->SolverBase<dim, degree, number>::solve();

  // If the FieldSolveType is constant or the solver is empty we can just return early.
  if (this->get_field_solve_type() == FieldSolveType::ExplicitConstant ||
      this->solver_is_empty())
    {
      return;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentSolver<dim, degree, number>::print()
{
  // Print the base class information
  this->SolverBase<dim, degree, number>::print();
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentSolver<dim, degree, number>::solve_explicit_equations(
  const std::function<
    void(std::vector<typename SolverBase<dim, degree, number>::VectorType *> &,
         const std::vector<typename SolverBase<dim, degree, number>::VectorType *> &)>
    &function)
{
  // Zero out the ghosts
  Timer::start_section("Zero ghosts");
  this->get_solution_handler().zero_out_ghosts();
  Timer::end_section("Zero ghosts");

  // Compute the update with the provided function
  function(new_solution_subset, solution_subset);

  // Scale the update by the respective (Scalar/Vector) invm. Note that we do this with
  // the original solution set to avoid some messy mapping.
  for (auto [index, vector] : this->get_solution_handler().get_new_solution_vector())
    {
      if (this->get_subset_attributes().find(index) !=
          this->get_subset_attributes().end())
        {
          vector->scale(this->get_invm_handler().get_invm(index));
        }
    }

  // Update the solutions
  this->get_solution_handler().update(this->get_field_solve_type(),
                                      this->get_solve_block());

  // Apply constraints
  // TODO (landinjm): This applies the constraints even to the old fields, which is
  // incorrect.
  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      this->get_solution_handler()
        .apply_constraints(index, this->get_constraint_handler().get_constraint(index));
    }

  // Update the ghosts
  Timer::start_section("Update ghosts");
  this->get_solution_handler().update_ghosts();
  Timer::end_section("Update ghosts");
}

#include "solvers/concurrent_solver.inst"

PRISMS_PF_END_NAMESPACE
