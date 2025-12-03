// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/timer.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/sequential_solver.h>
#include <prismspf/solvers/solver_base.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <map>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
SequentialSolver<dim, degree, number>::SequentialSolver(
  const SolverContext<dim, degree, number> &_solver_context,
  const FieldSolveType                     &_field_solve_type,
  Types::Index                              _solve_priority)
  : SolverBase<dim, degree, number>(_solver_context, _field_solve_type, _solve_priority)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSolver<dim, degree, number>::init()
{
  // Call the base class init
  this->SolverBase<dim, degree, number>::init();

  // If the solver is empty we can just return early.
  if (this->solver_is_empty())
    {
      return;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSolver<dim, degree, number>::reinit()
{
  // Call the base class reinit
  this->SolverBase<dim, degree, number>::reinit();

  // If the solver is empty we can just return early.
  if (this->solver_is_empty())
    {
      return;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSolver<dim, degree, number>::solve()
{
  // Call the base class solve
  this->SolverBase<dim, degree, number>::solve();

  // If the solver is empty we can just return early.
  if (this->solver_is_empty())
    {
      return;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSolver<dim, degree, number>::print()
{
  // Print the base class information
  this->SolverBase<dim, degree, number>::print();
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSolver<dim, degree, number>::init_linear_solver(
  const VariableAttributes &variable)
{
  // Grab the global field index
  const Types::Index global_field_index = variable.get_field_index();

  if (this->get_user_inputs()
        .get_linear_solve_parameters()
        .get_linear_solve_parameters(global_field_index)
        .preconditioner == PreconditionerType::GMG)
    {
      gmg_solvers.emplace(
        global_field_index,
        std::make_unique<GMGSolver<dim, degree, number>>(this->get_solver_context(),
                                                         variable));
      gmg_solvers.at(global_field_index)->init();
    }
  else
    {
      identity_solvers.emplace(
        global_field_index,
        std::make_unique<IdentitySolver<dim, degree, number>>(this->get_solver_context(),
                                                              variable));
      identity_solvers.at(global_field_index)->init();
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSolver<dim, degree, number>::init_explicit_solver(
  const VariableAttributes &variable)
{
  // Grab the global field index
  const Types::Index global_field_index = variable.get_field_index();

  // Creating temporary map to match types
  std::map<Types::Index, VariableAttributes> temp;
  temp.emplace(global_field_index, variable);
  subset_attributes_list.push_back(temp);

  // Create the implementation of MatrixFreeOperator with the subset of variable
  // attributes
  system_matrix[global_field_index] =
    std::make_unique<typename SolverBase<dim, degree, number>::SystemMatrixType>(
      subset_attributes_list.back(),
      this->get_pde_operator(),
      this->get_solve_block(),
      global_field_index);

  // Set up the user-implemented equations and create the residual vectors
  system_matrix[global_field_index]->clear();
  system_matrix[global_field_index]->initialize(
    this->get_matrix_free_container().get_matrix_free(),
    this->get_element_volume_container().get_element_volume());

  // Grab some data from the VariableAttributes
  const Types::Index max_fields =
    this->get_subset_attributes().begin()->second.get_max_fields();
  const Types::Index max_dependency_types =
    this->get_subset_attributes().begin()->second.get_max_dependency_types();

  // Resize the global to local solution vector
  global_to_local_solution[global_field_index].resize(
    static_cast<unsigned long>(max_fields) * max_dependency_types,
    Numbers::invalid_index);

  // Create the subset of solution vectors and add the mapping to MatrixFreeOperator
  new_solution_subset[global_field_index].push_back(
    this->get_solution_handler().get_new_solution_vector(global_field_index));
  solution_subset[global_field_index].push_back(
    this->get_solution_handler().get_solution_vector(global_field_index,
                                                     DependencyType::Normal));
  global_to_local_solution[global_field_index]
                          [(global_field_index * max_dependency_types) +
                           static_cast<Types::Index>(DependencyType::Normal)] = 0;

  Types::Index variable_index = 0;
  for (const auto &inner_dependency_set : variable.get_dependency_set_rhs())
    {
      Types::Index dependency_type = 0;
      for (const auto &field_type : inner_dependency_set)
        {
          // Skip if an invalid field type is found or the global_to_local_solution
          // already has an entry for this dependency index and dependency type
          if (field_type == FieldInfo::TensorRank::Undefined ||
              global_to_local_solution[global_field_index]
                                      [(variable_index * max_dependency_types) +
                                       dependency_type] != Numbers::invalid_index)
            {
              dependency_type++;
              continue;
            }

          solution_subset[global_field_index].push_back(
            this->get_solution_handler().get_solution_vector(variable_index,
                                                             static_cast<DependencyType>(
                                                               dependency_type)));
          global_to_local_solution[global_field_index][(variable_index *
                                                        max_dependency_types) +
                                                       dependency_type] =
            static_cast<unsigned int>(solution_subset.at(global_field_index).size()) - 1;

          dependency_type++;
        }

      variable_index++;
    }
  system_matrix[global_field_index]->add_global_to_local_mapping(
    global_to_local_solution[global_field_index]);
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSolver<dim, degree, number>::reinit_linear_solver(
  const VariableAttributes &variable)
{
  // Grab the global field index
  const Types::Index global_field_index = variable.get_field_index();

  if (this->get_user_inputs()
        .get_linear_solve_parameters()
        .get_linear_solve_parameters(global_field_index)
        .preconditioner == PreconditionerType::GMG)
    {
      gmg_solvers.at(global_field_index)->reinit();
    }
  else
    {
      identity_solvers.at(global_field_index)->reinit();
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSolver<dim, degree, number>::reinit_explicit_solver(
  [[maybe_unused]] const VariableAttributes &variable)
{
  // Grab the global field index
  const Types::Index global_field_index = variable.get_field_index();

  // Clear some objects
  global_to_local_solution[global_field_index].clear();
  solution_subset[global_field_index].clear();
  new_solution_subset[global_field_index].clear();

  // Set up the user-implemented equations and create the residual vectors
  system_matrix[global_field_index]->clear();
  system_matrix[global_field_index]->initialize(
    this->get_matrix_free_container().get_matrix_free(),
    this->get_element_volume_container().get_element_volume());

  // Grab some data from the VariableAttributes
  const Types::Index max_fields =
    this->get_subset_attributes().begin()->second.get_max_fields();
  const Types::Index max_dependency_types =
    this->get_subset_attributes().begin()->second.get_max_dependency_types();

  // Resize the global to local solution vector
  global_to_local_solution[global_field_index].resize(
    static_cast<unsigned long>(max_fields) * max_dependency_types,
    Numbers::invalid_index);

  // Create the subset of solution vectors and add the mapping to MatrixFreeOperator
  new_solution_subset[global_field_index].push_back(
    this->get_solution_handler().get_new_solution_vector(global_field_index));
  solution_subset[global_field_index].push_back(
    this->get_solution_handler().get_solution_vector(global_field_index,
                                                     DependencyType::Normal));
  global_to_local_solution[global_field_index]
                          [(global_field_index * max_dependency_types) +
                           static_cast<Types::Index>(DependencyType::Normal)] = 0;

  Types::Index variable_index = 0;
  for (const auto &inner_dependency_set : variable.get_dependency_set_rhs())
    {
      Types::Index dependency_type = 0;
      for (const auto &field_type : inner_dependency_set)
        {
          // Skip if an invalid field type is found or the global_to_local_solution
          // already has an entry for this dependency index and dependency type
          if (field_type == FieldInfo::TensorRank::Undefined ||
              global_to_local_solution[global_field_index]
                                      [(variable_index * max_dependency_types) +
                                       dependency_type] != Numbers::invalid_index)
            {
              dependency_type++;
              continue;
            }

          solution_subset[global_field_index].push_back(
            this->get_solution_handler().get_solution_vector(variable_index,
                                                             static_cast<DependencyType>(
                                                               dependency_type)));
          global_to_local_solution[global_field_index][(variable_index *
                                                        max_dependency_types) +
                                                       dependency_type] =
            static_cast<unsigned int>(solution_subset.at(global_field_index).size()) - 1;

          dependency_type++;
        }

      variable_index++;
    }
  system_matrix[global_field_index]->add_global_to_local_mapping(
    global_to_local_solution[global_field_index]);
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSolver<dim, degree, number>::solve_explicit_solver(
  const VariableAttributes &variable)
{
  // Zero out the ghosts
  Timer::start_section("Zero ghosts");
  this->get_solution_handler().zero_out_ghosts();
  Timer::end_section("Zero ghosts");

  // Grab the global field index
  const Types::Index global_field_index = variable.get_field_index();

  // Compute the update
  system_matrix[global_field_index]->compute_nonexplicit_auxiliary_update(
    new_solution_subset.at(global_field_index),
    solution_subset.at(global_field_index));

  // Scale the update by the respective (Scalar/Vector) invm.
  new_solution_subset.at(global_field_index)
    .at(0)
    ->scale(this->get_invm_handler().get_invm(global_field_index));

  // Update the solutions
  this->get_solution_handler().update(this->get_field_solve_type(),
                                      this->get_solve_block(),
                                      global_field_index);

  // Apply constraints
  this->get_constraint_handler()
    .get_constraint(global_field_index)
    .distribute(
      *(this->get_solution_handler().get_solution_vector(global_field_index,
                                                         DependencyType::Normal)));

  // Update the ghosts
  Timer::start_section("Update ghosts");
  this->get_solution_handler().update_ghosts();
  Timer::end_section("Update ghosts");
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialSolver<dim, degree, number>::solve_linear_solver(
  const VariableAttributes &variable)
{
  // Grab the old solution in a temporary variable
  const auto old_solution =
    *(this->get_solution_handler().get_solution_vector(variable.get_field_index(),
                                                       DependencyType::Normal));

  // Zero out the ghosts
  Timer::start_section("Zero ghosts");
  this->get_solution_handler().zero_out_ghosts();
  Timer::end_section("Zero ghosts");

  // Grab the global field index
  const Types::Index global_field_index = variable.get_field_index();

  // Skip if the field type is ImplicitTimeDependent and the increment is 0.
  if (variable.get_pde_type() == PDEType::ImplicitTimeDependent &&
      this->get_user_inputs().get_temporal_discretization().get_increment() == 0)
    {
      return;
    }

  if (this->get_user_inputs()
        .get_linear_solve_parameters()
        .get_linear_solve_parameters(global_field_index)
        .preconditioner == PreconditionerType::GMG)
    {
      gmg_solvers.at(global_field_index)->solve();
    }
  else
    {
      identity_solvers.at(global_field_index)->solve();
    }

  // The solve will have updated the "old solution" vector with the newton update so it's
  // technically the new solution. In order to update the solutions and perserve the old
  // states we copy the old solution from above and swap.
  *(this->get_solution_handler().get_new_solution_vector(variable.get_field_index())) =
    old_solution;
  this->get_solution_handler()
    .get_solution_vector(variable.get_field_index(), DependencyType::Normal)
    ->swap(
      *this->get_solution_handler().get_new_solution_vector(variable.get_field_index()));

  // Update the solutions
  this->get_solution_handler().update(this->get_field_solve_type(),
                                      this->get_solve_block(),
                                      global_field_index);

  // Update the ghosts
  Timer::start_section("Update ghosts");
  this->get_solution_handler().update_ghosts();
  Timer::end_section("Update ghosts");
}

template <unsigned int dim, unsigned int degree, typename number>
number
SequentialSolver<dim, degree, number>::solve_linear_solver(
  const VariableAttributes &variable,
  const number             &step_length)
{
  // Grab the global field index
  const Types::Index global_field_index = variable.get_field_index();

  // Skip if the field type is ImplicitTimeDependent and the increment is 0.
  if (variable.get_pde_type() == PDEType::ImplicitTimeDependent &&
      this->get_user_inputs().get_temporal_discretization().get_increment() == 0)
    {
      return 0.0;
    }

  if (this->get_user_inputs()
        .get_linear_solve_parameters()
        .get_linear_solve_parameters(global_field_index)
        .preconditioner == PreconditionerType::GMG)
    {
      gmg_solvers.at(global_field_index)->solve(step_length);
    }
  else
    {
      identity_solvers.at(global_field_index)->solve(step_length);
    }

  // Update the ghosts
  Timer::start_section("Update ghosts");
  this->get_solution_handler().update_ghosts();
  Timer::end_section("Update ghosts");

  // Return the norm of the newton update
  if (this->get_user_inputs()
        .get_linear_solve_parameters()
        .get_linear_solve_parameters(global_field_index)
        .preconditioner == PreconditionerType::GMG)
    {
      return gmg_solvers.at(global_field_index)->get_newton_update_l2_norm();
    }
  return identity_solvers.at(global_field_index)->get_newton_update_l2_norm();
}

#include "solvers/sequential_solver.inst"

PRISMS_PF_END_NAMESPACE
