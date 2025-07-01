// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/explicit_base.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <ostream>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
ExplicitBase<dim, degree>::ExplicitBase(const SolverContext<dim, degree> &_solver_context)
  : solver_context(&_solver_context)
{}

template <unsigned int dim, unsigned int degree>
void
ExplicitBase<dim, degree>::compute_subset_attributes(
  const FieldSolveType &field_solve_type)
{
  Assert((field_solve_type == FieldSolveType::Explicit ||
          field_solve_type == FieldSolveType::ExplicitPostprocess ||
          field_solve_type == FieldSolveType::ExplicitConstant),
         dealii::ExcMessage(
           "compute_subset_attributes() should only be used for "
           "Explicit, ExplicitPostprocess, and ExplicitConstant fieldSolveTypes"));

  subset_attributes.clear();

  for (const auto &[index, variable] :
       solver_context->get_user_inputs().get_variable_attributes())
    {
      if (variable.get_field_solve_type() == field_solve_type)
        {
          subset_attributes.emplace(index, variable);
        }
    }
}

template <unsigned int dim, unsigned int degree>
void
ExplicitBase<dim, degree>::compute_shared_dependencies()
{
  // Compute the shared dependency flags
  auto &dependency_flag_set = subset_attributes.begin()->second.get_eval_flag_set_rhs();
  for (const auto &[index, variable] : subset_attributes)
    {
      if (!variable.get_eval_flag_set_rhs().empty())
        {
          for (const auto &[pair, flag] : variable.get_eval_flag_set_rhs())
            {
              dependency_flag_set[pair] |= flag;
            }
        }
    }
  for (auto &[index, variable] : subset_attributes)
    {
      for (const auto &[pair, flag] : dependency_flag_set)
        {
          variable.get_eval_flag_set_rhs()[pair] |= flag;
        }
    }

  // Compute the shared dependency set
  auto &dependency_set = subset_attributes.begin()->second.get_dependency_set_rhs();
  for (const auto &[main_index, variable] : subset_attributes)
    {
      for (const auto &[dependency_index, map] : variable.get_dependency_set_rhs())
        {
          for (const auto &[dependency_type, field_type] : map)
            {
              dependency_set[dependency_index].emplace(dependency_type, field_type);
            }
        }
    }
  for (auto &[index, variable] : subset_attributes)
    {
      variable.set_dependency_set_rhs(dependency_set);
    }

#ifdef DEBUG
  print();
#endif
}

template <unsigned int dim, unsigned int degree>
void
ExplicitBase<dim, degree>::set_initial_condition()
{
  for (const auto &[index, variable] : subset_attributes)
    {
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
                    ReadInitialCondition<dim>(
                      initial_condition_file.filename + "." +
                        initial_condition_file.file_extension,
                      initial_condition_file.file_variable_names
                        [iterator -
                         initial_condition_file.simulation_variable_names.begin()],
                      subset_attributes.at(index).get_field_type()),
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
            InitialCondition<dim, degree>(index,
                                          subset_attributes.at(index).get_field_type(),
                                          solver_context->get_pde_operator()),
            *(solver_context->get_solution_handler()
                .get_solution_vector(index, DependencyType::Normal)));
        }
    }
}

template <unsigned int dim, unsigned int degree>
void
ExplicitBase<dim, degree>::print()
{
  ConditionalOStreams::pout_summary()
    << "  ==============================================\n"
    << "    Shared dependency set\n"
    << "  ==============================================\n";
  const auto &dependency_set = subset_attributes.begin()->second.get_dependency_set_rhs();
  for (const auto &[index, map] : dependency_set)
    {
      for (const auto &[dependency_type, field_type] : map)
        {
          ConditionalOStreams::pout_summary()
            << "  Index: " << index << " Dependency: " << to_string(dependency_type)
            << " Field: " << to_string(field_type) << "\n";
        }
    }
  ConditionalOStreams::pout_summary() << "\n" << std::flush;
}

INSTANTIATE_BI_TEMPLATE(ExplicitBase)

PRISMS_PF_END_NAMESPACE
