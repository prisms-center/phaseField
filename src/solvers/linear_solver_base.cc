// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/linear_solver_base.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <memory>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
LinearSolverBase<dim, degree, number>::LinearSolverBase(
  const SolverContext<dim, degree, number> &_solver_context,
  const VariableAttributes                 &_variable_attributes)
  : solver_context(&_solver_context)
  , variable_attributes(&_variable_attributes)
  , field_index(_variable_attributes.get_field_index())
  , residual(_solver_context.get_solution_handler().get_new_solution_vector(field_index))
  , newton_update(
      _solver_context.get_solution_handler().get_solution_vector(field_index,
                                                                 DependencyType::Change))
  , solver_control(_solver_context.get_user_inputs()
                     .get_linear_solve_parameters()
                     .get_linear_solve_parameters(field_index)
                     .max_iterations)
{
  // Creating map to match types
  subset_attributes.emplace(field_index, *variable_attributes);

  // Create the implementation of MatrixFreeOperator with the subset of variable
  // attributes
  system_matrix =
    std::make_unique<SystemMatrixType>(subset_attributes,
                                       _solver_context.get_pde_operator(),
                                       variable_attributes->get_solve_block(),
                                       field_index);
  update_system_matrix =
    std::make_unique<SystemMatrixType>(subset_attributes,
                                       _solver_context.get_pde_operator(),
                                       variable_attributes->get_solve_block(),
                                       field_index);

  // Grab some data from the VariableAttributes
  const Types::Index max_fields = variable_attributes->get_max_fields();
  const Types::Index max_dependency_types =
    variable_attributes->get_max_dependency_types();

  // Resize the global to local solution vector
  residual_global_to_local_solution.resize(max_fields * max_dependency_types,
                                           Numbers::invalid_index);

  newton_update_global_to_local_solution.resize(max_fields * max_dependency_types,
                                                Numbers::invalid_index);

  // Create the residual subset of solution vectors and add the mapping to
  // MatrixFreeOperator
  residual_src.push_back(
    _solver_context.get_solution_handler().get_solution_vector(field_index,
                                                               DependencyType::Normal));
  residual_global_to_local_solution[(field_index * max_dependency_types) +
                                    static_cast<Types::Index>(DependencyType::Normal)] =
    0;

  Types::Index variable_index = 0;
  for (const auto &inner_dependency_set : variable_attributes->get_dependency_set_rhs())
    {
      Types::Index dependency_type = 0;
      for (const auto &field_type : inner_dependency_set)
        {
          // Skip if an invalid field type is found or the global_to_local_solution
          // already has an entry for this dependency index and dependency type
          if (field_type == Numbers::invalid_field_type ||
              residual_global_to_local_solution[(variable_index * max_dependency_types) +
                                                dependency_type] !=
                Numbers::invalid_index)
            {
              dependency_type++;
              continue;
            }

          residual_src.push_back(
            _solver_context.get_solution_handler().get_solution_vector(
              variable_index,
              static_cast<DependencyType>(dependency_type)));
          residual_global_to_local_solution[(variable_index * max_dependency_types) +
                                            dependency_type] = residual_src.size() - 1;

          dependency_type++;
        }

      variable_index++;
    }

  // Create the newton update subset of solution vectors and add the mapping to
  // MatrixFreeOperator. For this one we consider the singular src vector as the residual
  // vector above and the src subset all other dependencies vectors. This is
  // complicated, but has to do with the way deal.II's solvers (like CG) handle
  // iterative updates. For this reason, we have to pass the residual vector as
  // VectorType src and all other dependencies for the LHS as std::vector<VectorType*>
  // src_subset.

  variable_index = 0;
  for (const auto &inner_dependency_set : variable_attributes->get_dependency_set_lhs())
    {
      Types::Index dependency_type = 0;
      for (const auto &field_type : inner_dependency_set)
        {
          // Skip if an invalid field type is found or the global_to_local_solution
          // already has an entry for this dependency index and dependency type
          if (field_type == Numbers::invalid_field_type ||
              newton_update_global_to_local_solution
                  [(variable_index * max_dependency_types) + dependency_type] !=
                Numbers::invalid_index)
            {
              dependency_type++;
              continue;
            }

          if (static_cast<DependencyType>(dependency_type) == DependencyType::Change)
            {
              Assert(field_index == variable_index,
                     dealii::ExcMessage("The change type should have the same type as "
                                        "the field we're solving."));
            }
          newton_update_src.push_back(
            _solver_context.get_solution_handler().get_solution_vector(
              variable_index,
              static_cast<DependencyType>(dependency_type)));
          newton_update_global_to_local_solution[(variable_index * max_dependency_types) +
                                                 dependency_type] =
            newton_update_src.size() - 1;

          dependency_type++;
        }

      variable_index++;
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
LinearSolverBase<dim, degree, number>::compute_solver_tolerance()
{
  tolerance = solver_context->get_user_inputs()
                    .get_linear_solve_parameters()
                    .get_linear_solve_parameters(field_index)
                    .tolerance_type == SolverToleranceType::RelativeResidualChange
                ? solver_context->get_user_inputs()
                      .get_linear_solve_parameters()
                      .get_linear_solve_parameters(field_index)
                      .tolerance *
                    residual->l2_norm()
                : solver_context->get_user_inputs()
                    .get_linear_solve_parameters()
                    .get_linear_solve_parameters(field_index)
                    .tolerance;
}

#include "solvers/linear_solver_base.inst"

PRISMS_PF_END_NAMESPACE
