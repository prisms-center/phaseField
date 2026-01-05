// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/conditional_ostreams.h>

#include <prismspf/solvers/group_solver_base.h>

#include "prismspf/core/group_solution_handler.h"

#include <map>
#include <ostream>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
GroupSolverBase<dim, degree, number>::GroupSolverBase(
  SolveGroup                                _solve_group,
  const SolverContext<dim, degree, number> &_solver_context)
  : solve_group(std::move(_solve_group))
  , solutions()
  , solver_context(std::make_shared<SolverContext<dim, degree, number>>(_solver_context))
{}

template <unsigned int dim, unsigned int degree, typename number>
void
GroupSolverBase<dim, degree, number>::init()
{
  // Initialize vectors
  solutions =
    GroupSolutionHandler<dim, number>(solve_group, solver_context->field_attributes);
  solutions.init(solver_context->mapping,
                 solver_context->dof_manager,
                 solver_context->constraint_manager,
                 solver_context->quadrature);

  // Set the initial condition
  set_initial_condition();

  // Apply constraints.
  solutions.apply_constraints();

  // Initialize mf_operators
}

template <unsigned int dim, unsigned int degree, typename number>
void
GroupSolverBase<dim, degree, number>::reinit()
{
  // Apply constraints.
  solutions.reinit();
  solutions.apply_constraints();
}

template <unsigned int dim, unsigned int degree, typename number>
void
GroupSolverBase<dim, degree, number>::solve()
{}

template <unsigned int dim, unsigned int degree, typename number>
void
GroupSolverBase<dim, degree, number>::solve_level(
  [[maybe_unused]] unsigned int relative_level)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
GroupSolverBase<dim, degree, number>::print()
{}

template <unsigned int dim, unsigned int degree, typename number>
void
GroupSolverBase<dim, degree, number>::set_initial_condition()
{
  for (const auto &global_index : solve_group.field_indices)
    {
      if (solver_context->get_user_inputs()
            .get_load_initial_condition_parameters()
            .get_read_initial_conditions_from_file())
        {
          const auto &initial_condition_parameters =
            solver_context->get_user_inputs().get_load_initial_condition_parameters();
          for (const auto &initial_condition_file :
               initial_condition_parameters.get_initial_condition_files())
            {
              auto name_it =
                std::find(initial_condition_file.simulation_variable_names.begin(),
                          initial_condition_file.simulation_variable_names.end(),
                          attributes_list[global_index].name);
              if (name_it != initial_condition_file.simulation_variable_names.end())
                {
                  dealii::VectorTools::interpolate(
                    solver_context->get_mapping(),
                    solver_context->get_dof_manager().get_dof_handler(global_index),
                    ReadInitialCondition<dim, number>(
                      *name_it,
                      attributes_list[global_index].field_type,
                      initial_condition_file,
                      solver_context->get_user_inputs().get_spatial_discretization()),
                    solutions.get_solution_vector(global_index));
                }
            }
        }
      else
        {
          dealii::VectorTools::interpolate(
            solver_context->get_mapping(),
            solver_context->get_dof_manager().get_dof_handler(index),
            InitialCondition<dim, degree, number>(
              index,
              attributes_list[global_index].field_type,
              solver_context->get_pde_operator()),
            solutions.get_solution_vector(global_index));
        }
      solutions.apply_initial_condition_for_old_fields();
    }
}

// #include "solvers/group_solver_base.inst"

PRISMS_PF_END_NAMESPACE