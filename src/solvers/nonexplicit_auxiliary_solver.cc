// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/timer.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/solvers/nonexplicit_auxiliary_solver.h>
#include <prismspf/solvers/nonexplicit_base.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <map>
#include <memory>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
NonexplicitAuxiliarySolver<dim, degree>::NonexplicitAuxiliarySolver(
  const SolverContext<dim, degree> &_solver_context)
  : NonexplicitBase<dim, degree>(_solver_context)
{}

template <unsigned int dim, unsigned int degree>
void
NonexplicitAuxiliarySolver<dim, degree>::init()
{
  this->compute_subset_attributes(FieldSolveType::NonexplicitAuxiliary);

  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }

  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      // Creating temporary map to match types
      std::map<unsigned int, VariableAttributes> temp;
      temp.emplace(index, variable);
      subset_attributes_list.push_back(temp);

      // Create the implementation of MatrixFreeOperator with the subset of variable
      // attributes
      this->get_system_matrix()[index] =
        std::make_unique<SystemMatrixType>(this->get_subset_attributes(),
                                           this->get_pde_operator(),
                                           index);

      // Set up the user-implemented equations and create the residual vectors
      this->get_system_matrix()[index]->clear();
      this->get_system_matrix()[index]->initialize(
        this->get_matrix_free_handler().get_matrix_free());

      // Resize the global to local solution vector
      global_to_local_solution[index].resize(variable.get_dependency_set_rhs().size());
      for (auto &vector : global_to_local_solution[index])
        {
          vector.resize(variable.get_dependency_set_rhs().begin()->size(),
                        Numbers::invalid_index);
        }

      // Create the subset of solution vectors and add the mapping to MatrixFreeOperator
      new_solution_subset[index].push_back(
        this->get_solution_handler().get_new_solution_vector(index));
      solution_subset[index].push_back(
        this->get_solution_handler().get_solution_vector(index, DependencyType::Normal));
      global_to_local_solution[index][index]
                              [static_cast<Types::Index>(DependencyType::Normal)] = 0;

      Types::Index variable_index = 0;
      for (const auto &inner_dependency_set :
           this->get_subset_attributes().begin()->second.get_dependency_set_rhs())
        {
          Types::Index dependency_type = 0;
          for (const auto &field_type : inner_dependency_set)
            {
              // Skip if an invalid field type is found or the global_to_local_solution
              // already has an entry for this dependency index and dependency type
              if (field_type == Numbers::invalid_field_type ||
                  global_to_local_solution[index][variable_index][dependency_type] !=
                    Numbers::invalid_index)
                {
                  dependency_type++;
                  continue;
                }

              solution_subset[index].push_back(
                this->get_solution_handler().get_solution_vector(
                  variable_index,
                  static_cast<DependencyType>(dependency_type)));
              global_to_local_solution[index][variable_index][dependency_type] =
                solution_subset.at(index).size() - 1;

              dependency_type++;
            }

          variable_index++;
        }
      this->get_system_matrix()[index]->add_global_to_local_mapping(
        global_to_local_solution[index]);
    }
}

template <unsigned int dim, unsigned int degree>
void
NonexplicitAuxiliarySolver<dim, degree>::solve()
{
  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }

  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      // Compute the update
      this->get_system_matrix()[index]->compute_nonexplicit_auxiliary_update(
        new_solution_subset.at(index),
        solution_subset.at(index));

      // Scale the update by the respective (Scalar/Vector) invm.
      new_solution_subset.at(index).at(0)->scale(
        this->get_invm_handler().get_invm(index));

      // Update the solutions
      this->get_solution_handler().update(FieldSolveType::NonexplicitAuxiliary, index);

      // Apply constraints
      this->get_constraint_handler().get_constraint(index).distribute(*(
        this->get_solution_handler().get_solution_vector(index, DependencyType::Normal)));

      // Update the ghosts
      Timer::start_section("Update ghosts");
      this->get_solution_handler().update_ghosts();
      Timer::end_section("Update ghosts");
    }
}

INSTANTIATE_BI_TEMPLATE(NonexplicitAuxiliarySolver)

PRISMS_PF_END_NAMESPACE