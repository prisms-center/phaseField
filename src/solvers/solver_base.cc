// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/solver_base.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <map>
#include <ostream>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
SolverBase<dim, degree, number>::SolverBase(
  const SolverContext<dim, degree, number> &_solver_context,
  const FieldSolveType                     &_field_solve_type,
  Types::Index                              _solve_priority)
  : solver_context(std::make_shared<SolverContext<dim, degree, number>>(_solver_context))
  , field_solve_type(_field_solve_type)
  , solve_priority(_solve_priority)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverBase<dim, degree, number>::init()
{
  // Update the subset of variable attributes
  update_subset_attributes(field_solve_type, solve_priority);

  // If the subset attribute is empty return early
  if (solver_is_empty())
    {
      ConditionalOStreams::pout_base() << "  no fields for this solver exist\n"
                                       << std::flush;
      return;
    }

  // Set the initial condition
  set_initial_condition();

  // Apply constraints. This part is neccessary so they are taken into account for
  // adaptive meshing
  for (const auto &[index, variable] : subset_attributes)
    {
      get_constraint_handler().get_constraint(index).distribute(
        *(get_solution_handler().get_solution_vector(index, DependencyType::Normal)));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverBase<dim, degree, number>::reinit()
{
  // If the subset attribute is empty return early
  if (solver_is_empty())
    {
      return;
    }

  // Apply constraints. This part is neccessary so they are taken into account for
  // adaptive meshing
  for (const auto &[index, variable] : subset_attributes)
    {
      get_constraint_handler().get_constraint(index).distribute(
        *(get_solution_handler().get_solution_vector(index, DependencyType::Normal)));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverBase<dim, degree, number>::solve()
{
  // If the subset attribute is empty return early
  if (solver_is_empty())
    {
      return;
    }

  // Do nothing
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverBase<dim, degree, number>::print()
{}

template <unsigned int dim, unsigned int degree, typename number>
bool
SolverBase<dim, degree, number>::solver_is_empty() const
{
  return subset_attributes.empty();
}

template <unsigned int dim, unsigned int degree, typename number>
std::map<Types::Index, VariableAttributes>
SolverBase<dim, degree, number>::compute_subset_attributes(
  const FieldSolveType &_field_solve_type,
  Types::Index          _solve_priority) const
{
  std::map<Types::Index, VariableAttributes> local_subset_attributes;

  for (const auto &[index, variable] :
       solver_context->get_user_inputs().get_variable_attributes())
    {
      if (variable.get_field_solve_type() == _field_solve_type &&
          variable.get_solve_block() == _solve_priority)
        {
          local_subset_attributes.emplace(index, variable);
        }
    }

  return local_subset_attributes;
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverBase<dim, degree, number>::update_subset_attributes(
  const FieldSolveType &_field_solve_type,
  Types::Index          _solve_priority)
{
  subset_attributes = compute_subset_attributes(_field_solve_type, _solve_priority);
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverBase<dim, degree, number>::set_initial_condition()
{
  for (const auto &[index, variable] : subset_attributes)
    {
      // TODO (landinjm): Skip certain fields for initial conditions

      Assert(solver_context->get_dof_handler().get_dof_handlers().size() > index,
             dealii::ExcMessage(
               "The const DoFHandler set is smaller than the given index = " +
               std::to_string(index)));
      Assert(subset_attributes.contains(index),
             dealii::ExcMessage(
               "There is no entry in the attribute subset for the given index = " +
               std::to_string(index)));

      if (solver_context->get_user_inputs()
            .get_load_initial_condition_parameters()
            .get_read_initial_conditions_from_file())
        {
          auto &initial_condition_parameters =
            solver_context->get_user_inputs().get_load_initial_condition_parameters();
          for (const auto &initial_condition_file :
               initial_condition_parameters.get_initial_condition_files())
            {
              auto iterator =
                std::find(initial_condition_file.simulation_variable_names.begin(),
                          initial_condition_file.simulation_variable_names.end(),
                          variable.get_name());
              if (iterator != initial_condition_file.simulation_variable_names.end())
                {
                  dealii::VectorTools::interpolate(
                    solver_context->get_mapping(),
                    *(solver_context->get_dof_handler().get_dof_handlers().at(index)),
                    ReadInitialCondition<dim, number>(
                      initial_condition_file.file_variable_names
                        [iterator -
                         initial_condition_file.simulation_variable_names.begin()],
                      subset_attributes.at(index).field_info.tensor_rank,
                      initial_condition_file,
                      solver_context->get_user_inputs().get_spatial_discretization()),
                    *(solver_context->get_solution_handler()
                        .get_solution_vector(index, DependencyType::Normal)));
                }
            }
        }
      else
        {
          dealii::VectorTools::interpolate(
            solver_context->get_mapping(),
            *(solver_context->get_dof_handler().get_dof_handlers().at(index)),
            InitialCondition<dim, degree, number>(
              index,
              subset_attributes.at(index).field_info.tensor_rank,
              solver_context->get_pde_operator()),
            *(solver_context->get_solution_handler()
                .get_solution_vector(index, DependencyType::Normal)));
        }

      // TODO (landinjm): Fix so that we apply some sort of initial condition to all old
      // vector for all types.
      solver_context->get_solution_handler().apply_initial_condition_for_old_fields();
    }
}

#include "solvers/solver_base.inst"

PRISMS_PF_END_NAMESPACE