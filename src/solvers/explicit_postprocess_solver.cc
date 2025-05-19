// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/explicit_base.h>
#include <prismspf/solvers/explicit_postprocess_solver.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
explicitPostprocessSolver<dim, degree>::explicitPostprocessSolver(
  const UserInputParameters<dim>                         &_user_inputs,
  const MatrixfreeHandler<dim>                           &_matrix_free_handler,
  const InvmHandler<dim, degree>                         &_invm_handler,
  const ConstraintHandler<dim, degree>                   &_constraint_handler,
  const DofHandler<dim>                                  &_dof_handler,
  const dealii::MappingQ1<dim>                           &_mapping,
  SolutionHandler<dim>                                   &_solution_handler,
  std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator)
  : explicitBase<dim, degree>(_user_inputs,
                              _matrix_free_handler,
                              _invm_handler,
                              _constraint_handler,
                              _dof_handler,
                              _mapping,
                              _solution_handler,
                              std::move(_pde_operator))
{}

template <unsigned int dim, unsigned int degree>
void
explicitPostprocessSolver<dim, degree>::init()
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

  // Create the subset of solution vectors and add the mapping to MatrixFreeOperator
  for (const auto &[index, map] :
       this->get_subset_attributes().begin()->second.get_dependency_set_rhs())
    {
      for (const auto &[dependency_type, field_type] : map)
        {
          const auto pair = std::make_pair(index, dependency_type);

          solution_subset.push_back(
            this->get_solution_handler().get_solution_vector(index, dependency_type));
          new_solution_subset.push_back(
            this->get_solution_handler().get_new_solution_vector(index));
          global_to_local_solution.emplace(pair, solution_subset.size() - 1);
        }
    }

  this->get_system_matrix()->add_global_to_local_mapping(global_to_local_solution);
}

template <unsigned int dim, unsigned int degree>
void
explicitPostprocessSolver<dim, degree>::solve()
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
}

INSTANTIATE_BI_TEMPLATE(explicitPostprocessSolver)

PRISMS_PF_END_NAMESPACE