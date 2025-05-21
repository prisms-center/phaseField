// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/nonexplicit_base.h>

#include <prismspf/config.h>

#include <memory>
#include <ostream>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
NonexplicitBase<dim, degree>::NonexplicitBase(
  const UserInputParameters<dim>                         &_user_inputs,
  const MatrixfreeHandler<dim>                           &_matrix_free_handler,
  const TriangulationHandler<dim>                        &_triangulation_handler,
  const InvmHandler<dim, degree>                         &_invm_handler,
  const ConstraintHandler<dim, degree>                   &_constraint_handler,
  const DofHandler<dim>                                  &_dof_handler,
  const dealii::MappingQ1<dim>                           &_mapping,
  dealii::MGLevelObject<MatrixfreeHandler<dim, float>>   &_mg_matrix_free_handler,
  SolutionHandler<dim>                                   &_solution_handler,
  std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator)
  : user_inputs(&_user_inputs)
  , matrix_free_handler(&_matrix_free_handler)
  , triangulation_handler(&_triangulation_handler)
  , invm_handler(&_invm_handler)
  , constraint_handler(&_constraint_handler)
  , dof_handler(&_dof_handler)
  , mapping(&_mapping)
  , mg_matrix_free_handler(&_mg_matrix_free_handler)
  , solution_handler(&_solution_handler)
  , pde_operator(std::move(_pde_operator))
{}

template <unsigned int dim, unsigned int degree>
void
NonexplicitBase<dim, degree>::compute_subset_attributes(
  const FieldSolveType &field_solve_type)
{
  Assert((field_solve_type == FieldSolveType::NonexplicitLinear ||
          field_solve_type == FieldSolveType::NonexplicitSelfnonlinear ||
          field_solve_type == FieldSolveType::NonexplicitAuxiliary ||
          field_solve_type == FieldSolveType::NonexplicitCononlinear),
         dealii::ExcMessage(
           "compute_subset_attributes() should only be used for "
           "NonexplicitLinear, NonexplicitSelfnonlinear, NonexplicitAuxiliary, and "
           "NonexplicitCononlinear fieldSolveTypes"));

  subset_attributes.clear();

  for (const auto &[index, variable] : user_inputs->get_variable_attributes())
    {
      if (variable.get_field_solve_type() == field_solve_type)
        {
          subset_attributes.emplace(index, variable);
        }
    }
}

template <unsigned int dim, unsigned int degree>
void
NonexplicitBase<dim, degree>::compute_shared_dependencies()
{
  Assert(subset_attributes.begin()->second.get_field_solve_type() ==
           FieldSolveType::NonexplicitCononlinear,
         dealii::ExcMessage("compute_shared_dependencies() should only be used for "
                            "NonexplicitCononlinear fieldSolveTypes"));

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
NonexplicitBase<dim, degree>::set_initial_condition()
{
  for (const auto &[index, variable] : subset_attributes)
    {
      if (variable.get_pde_type() != PDEType::ImplicitTimeDependent &&
          variable.get_pde_type() != PDEType::TimeIndependent)
        {
          continue;
        }

      Assert(dof_handler->get_dof_handlers().size() > index,
             dealii::ExcMessage(
               "The const DoFHandler set is smaller than the given index = " +
               std::to_string(index)));
      Assert(subset_attributes.contains(index),
             dealii::ExcMessage(
               "There is no entry in the attribute subset for the given index = " +
               std::to_string(index)));

      dealii::VectorTools::interpolate(
        *mapping,
        *(dof_handler->get_dof_handlers().at(index)),
        InitialCondition<dim, degree>(index,
                                      subset_attributes.at(index).get_field_type(),
                                      pde_operator),
        *(solution_handler->get_solution_vector(index, DependencyType::Normal)));

      // TODO (landinjm): Fix so that we apply some sort of initial condition to all old
      // vector for all types.
      solution_handler->apply_initial_condition_for_old_fields();
    }
}

template <unsigned int dim, unsigned int degree>
void
NonexplicitBase<dim, degree>::print()
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

INSTANTIATE_BI_TEMPLATE(NonexplicitBase)

PRISMS_PF_END_NAMESPACE
