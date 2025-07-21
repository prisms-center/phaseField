// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/timer.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/solvers/explicit_base.h>
#include <prismspf/solvers/explicit_postprocess_solver.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
ExplicitPostprocessSolver<dim, degree>::ExplicitPostprocessSolver(
  const SolverContext<dim, degree> &_solver_context)
  : ExplicitBase<dim, degree>(_solver_context)
{}

template <unsigned int dim, unsigned int degree>
void
ExplicitPostprocessSolver<dim, degree>::init()
{
  this->compute_subset_attributes(FieldSolveType::ExplicitPostprocess);

  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }

  this->compute_shared_dependencies();

  // Create the implementation of MatrixFreeOperator with the subset of variable
  // attributes
  this->get_system_matrix() =
    std::make_unique<SystemMatrixType>(this->get_subset_attributes(),
                                       this->get_pde_operator());

  // Set up the user-implemented equations and create the residual vectors
  this->get_system_matrix()->clear();
  this->get_system_matrix()->initialize(
    this->get_matrix_free_handler().get_matrix_free());

  // Resize the global to local solution vector
  global_to_local_solution.resize(
    this->get_subset_attributes().begin()->second.get_dependency_set_rhs().size());
  for (auto &vector : global_to_local_solution)
    {
      vector.resize(this->get_subset_attributes()
                      .begin()
                      ->second.get_dependency_set_rhs()
                      .begin()
                      ->size(),
                    Numbers::invalid_index);
    }

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
          if (field_type == Numbers::invalid_field_type ||
              global_to_local_solution[dependency_index][dependency_type] !=
                Numbers::invalid_index)
            {
              dependency_type++;
              continue;
            }

          solution_subset.push_back(this->get_solution_handler().get_solution_vector(
            dependency_index,
            static_cast<DependencyType>(dependency_type)));
          new_solution_subset.push_back(
            this->get_solution_handler().get_new_solution_vector(dependency_index));
          global_to_local_solution[dependency_index][dependency_type] =
            solution_subset.size() - 1;

          dependency_type++;
        }

      dependency_index++;
    }

  this->get_system_matrix()->add_global_to_local_mapping(global_to_local_solution);
}

template <unsigned int dim, unsigned int degree>
void
ExplicitPostprocessSolver<dim, degree>::solve()
{
  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }

  // Compute the postprocessed fields
  this->get_system_matrix()->compute_postprocess_explicit_update(new_solution_subset,
                                                                 solution_subset);

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
  this->get_solution_handler().update(FieldSolveType::ExplicitPostprocess);

  // Update the ghosts
  Timer::start_section("Update ghosts");
  this->get_solution_handler().update_ghosts();
  Timer::end_section("Update ghosts");
}

INSTANTIATE_BI_TEMPLATE(ExplicitPostprocessSolver)

PRISMS_PF_END_NAMESPACE
