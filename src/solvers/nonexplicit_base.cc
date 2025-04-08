#include <prismspf/solvers/nonexplicit_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, int degree>
nonexplicitBase<dim, degree>::nonexplicitBase(
  const userInputParameters<dim>                         &_user_inputs,
  const matrixfreeHandler<dim>                           &_matrix_free_handler,
  const triangulationHandler<dim>                        &_triangulation_handler,
  const invmHandler<dim, degree>                         &_invm_handler,
  const constraintHandler<dim>                           &_constraint_handler,
  const dofHandler<dim>                                  &_dof_handler,
  const dealii::MappingQ1<dim>                           &_mapping,
  dealii::MGLevelObject<matrixfreeHandler<dim, float>>   &_mg_matrix_free_handler,
  solutionHandler<dim>                                   &_solution_handler,
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
  , pde_operator(_pde_operator)
{}

template <int dim, int degree>
inline void
nonexplicitBase<dim, degree>::compute_subset_attributes(
  const fieldSolveType &field_solve_type)
{
  Assert((field_solve_type == fieldSolveType::NONEXPLICIT_LINEAR ||
          field_solve_type == fieldSolveType::NONEXPLICIT_SELF_NONLINEAR ||
          field_solve_type == fieldSolveType::NONEXPLICIT_AUXILIARY ||
          field_solve_type == fieldSolveType::NONEXPLICIT_CO_NONLINEAR),
         dealii::ExcMessage(
           "compute_subset_attributes() should only be used for "
           "NONEXPLICIT_LINEAR, NONEXPLICIT_SELF_NONLINEAR, NONEXPLICIT_AUXILIARY, and "
           "NONEXPLICIT_CO_NONLINEAR fieldSolveTypes"));

  subset_attributes.clear();

  for (const auto &[index, variable] : *user_inputs->var_attributes)
    {
      if (variable.field_solve_type == field_solve_type)
        {
          subset_attributes.emplace(index, variable);
        }
    }
}

template <int dim, int degree>
inline void
nonexplicitBase<dim, degree>::compute_shared_dependencies()
{
  Assert(subset_attributes.begin()->second.field_solve_type ==
           fieldSolveType::NONEXPLICIT_CO_NONLINEAR,
         dealii::ExcMessage("compute_shared_dependencies() should only be used for "
                            "NONEXPLICIT_CO_NONLINEAR fieldSolveTypes"));

  // Compute the shared dependency flags
  auto &dependency_flag_set = subset_attributes.begin()->second.eval_flag_set_RHS;
  for (const auto &[index, variable] : subset_attributes)
    {
      if (!variable.eval_flag_set_RHS.empty())
        {
          for (const auto &[pair, flag] : variable.eval_flag_set_RHS)
            {
              dependency_flag_set[pair] |= flag;
            }
        }
    }
  for (auto &[index, variable] : subset_attributes)
    {
      for (const auto &[pair, flag] : dependency_flag_set)
        {
          variable.eval_flag_set_RHS[pair] |= flag;
        }
    }

  // Compute the shared dependency set
  auto &dependency_set = subset_attributes.begin()->second.dependency_set_RHS;
  for (const auto &[main_index, variable] : subset_attributes)
    {
      for (const auto &[dependency_index, map] : variable.dependency_set_RHS)
        {
          for (const auto &[dependency_type, field_type] : map)
            {
              dependency_set[dependency_index].emplace(dependency_type, field_type);
            }
        }
    }
  for (auto &[index, variable] : subset_attributes)
    {
      variable.dependency_set_RHS = dependency_set;
    }

#ifdef DEBUG
  print();
#endif
}

template <int dim, int degree>
inline void
nonexplicitBase<dim, degree>::set_initial_condition()
{
  for (const auto &[index, variable] : subset_attributes)
    {
      if (variable.pde_type != PDEType::IMPLICIT_TIME_DEPENDENT &&
          variable.pde_type != PDEType::TIME_INDEPENDENT)
        {
          continue;
        }

      Assert(dof_handler->get_dof_handlers().size() > index,
             dealii::ExcMessage(
               "The const DoFHandler set is smaller than the given index = " +
               std::to_string(index)));
      Assert(subset_attributes.find(index) != subset_attributes.end(),
             dealii::ExcMessage(
               "There is no entry in the attribute subset for the given index = " +
               std::to_string(index)));

      dealii::VectorTools::interpolate(
        *mapping,
        *(dof_handler->get_dof_handlers().at(index)),
        initialCondition<dim>(index,
                              subset_attributes.at(index).field_type,
                              *user_inputs),
        *(solution_handler->get_solution_vector(index, dependencyType::NORMAL)));

      // TODO (landinjm): Fix so that we apply some sort of initial condition to all old
      // vector for all types.
      solution_handler->apply_initial_condition_for_old_fields();
    }
}

template <int dim, int degree>
inline void
nonexplicitBase<dim, degree>::print()
{
  conditionalOStreams::pout_summary()
    << "  ==============================================\n"
    << "    Shared dependency set\n"
    << "  ==============================================\n";
  const auto &dependency_set = subset_attributes.begin()->second.dependency_set_RHS;
  for (const auto &[index, map] : dependency_set)
    {
      for (const auto &[dependency_type, field_type] : map)
        {
          conditionalOStreams::pout_summary()
            << "  Index: " << index << " Dependency: " << to_string(dependency_type)
            << " Field: " << to_string(field_type) << "\n";
        }
    }
  conditionalOStreams::pout_summary() << "\n" << std::flush;
}

INSTANTIATE_BI_TEMPLATE(nonexplicitBase)

PRISMS_PF_END_NAMESPACE