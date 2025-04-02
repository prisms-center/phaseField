// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/numerics/vector_tools.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * Forward declaration for user-implemented PDE class.
 */
template <int dim, int degree, typename number>
class customPDE;

/**
 * \brief Base class for nonexplicit solves.
 */
template <int dim, int degree>
class nonexplicitBase
{
public:
  using SystemMatrixType = customPDE<dim, degree, double>;

  /**
   * \brief Constructor.
   */
  nonexplicitBase(
    const userInputParameters<dim>                       &_user_inputs,
    const matrixfreeHandler<dim>                         &_matrix_free_handler,
    const triangulationHandler<dim>                      &_triangulation_handler,
    const invmHandler<dim, degree>                       &_invm_handler,
    const constraintHandler<dim>                         &_constraint_handler,
    const dofHandler<dim>                                &_dof_handler,
    const dealii::MappingQ1<dim>                         &_mapping,
    dealii::MGLevelObject<matrixfreeHandler<dim, float>> &_mg_matrix_free_handler,
    solutionHandler<dim>                                 &_solution_handler);

  /**
   * \brief Destructor.
   */
  ~nonexplicitBase() = default;

  /**
   * \brief Initialize system.
   */
  virtual void
  init() = 0;

  /**
   * \brief Solve a single update step.
   */
  virtual void
  solve() = 0;

protected:
  /**
   * \brief Compute the subset of variableAttributes that belongs to a given
   * fieldSolveType. This function should only be used for nonexplicit fieldSolveTypes,
   * such as NONEXPLICIT_LINEAR, NONEXPLICIT_SELF_NONLINEAR, NONEXPLICIT_AUXILIARY, and
   * NONEXPLICIT_CO_NONLINEAR.
   */
  void
  compute_subset_attributes(const fieldSolveType &field_solve_type);

  /**
   * \brief Compute the shared dependency set and copy it to all eval_flag_set_RHS. Also
   * do something similar with dependency_set_RHS so that all the FEEvaluation objects are
   * initialized. This should only be called for concurrent nonexplicit fieldSolveTypes
   * like NONEXPLICIT_CO_NONLINEAR.
   */
  void
  compute_shared_dependencies();

  /**
   * \brief Set the initial condition according to subset_attributes. This only applies
   * for PDEType IMPLICIT_TIME_DEPENDENT fields.
   */
  void
  set_initial_condition();

  /**
   * \brief Print dependency_set_RHS to summary.log
   */
  void
  print();

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> &user_inputs;

  /**
   * \brief Matrix-free object handler for non-multigrid data.
   */
  const matrixfreeHandler<dim> &matrix_free_handler;

  /**
   * \brief Triangulation handler.
   */
  const triangulationHandler<dim> &triangulation_handler;

  /**
   * \brief invm handler.
   */
  const invmHandler<dim, degree> &invm_handler;

  /**
   * \brief Constraint handler.
   */
  const constraintHandler<dim> &constraint_handler;

  /**
   * \brief DoF handler.
   */
  const dofHandler<dim> &dof_handler;

  /**
   * \brief Mappings to and from reference cell.
   */
  const dealii::MappingQ1<dim> &mapping;

  /**
   * \brief Matrix-free object handler for multigrid data.
   */
  dealii::MGLevelObject<matrixfreeHandler<dim, float>> &mg_matrix_free_handler;

  /**
   * \brief Solution handler.
   */
  solutionHandler<dim> &solution_handler;

  /**
   * \brief Subset of variable attributes for fields.
   */
  std::map<unsigned int, variableAttributes> subset_attributes;

  /**
   * \brief PDE operator for the residual side.
   */
  std::map<unsigned int, std::unique_ptr<SystemMatrixType>> system_matrix;

  /**
   * \brief PDE operator for the newton update side.
   */
  std::map<unsigned int, std::unique_ptr<SystemMatrixType>> update_system_matrix;
};

template <int dim, int degree>
nonexplicitBase<dim, degree>::nonexplicitBase(
  const userInputParameters<dim>                       &_user_inputs,
  const matrixfreeHandler<dim>                         &_matrix_free_handler,
  const triangulationHandler<dim>                      &_triangulation_handler,
  const invmHandler<dim, degree>                       &_invm_handler,
  const constraintHandler<dim>                         &_constraint_handler,
  const dofHandler<dim>                                &_dof_handler,
  const dealii::MappingQ1<dim>                         &_mapping,
  dealii::MGLevelObject<matrixfreeHandler<dim, float>> &_mg_matrix_free_handler,
  solutionHandler<dim>                                 &_solution_handler)
  : user_inputs(_user_inputs)
  , matrix_free_handler(_matrix_free_handler)
  , triangulation_handler(_triangulation_handler)
  , invm_handler(_invm_handler)
  , constraint_handler(_constraint_handler)
  , dof_handler(_dof_handler)
  , mapping(_mapping)
  , mg_matrix_free_handler(_mg_matrix_free_handler)
  , solution_handler(_solution_handler)
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

  for (const auto &[index, variable] : user_inputs.var_attributes)
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

      Assert(dof_handler.get_dof_handlers().size() > index,
             dealii::ExcMessage(
               "The const DoFHandler set is smaller than the given index = " +
               std::to_string(index)));
      Assert(subset_attributes.find(index) != subset_attributes.end(),
             dealii::ExcMessage(
               "There is no entry in the attribute subset for the given index = " +
               std::to_string(index)));

      dealii::VectorTools::interpolate(
        mapping,
        *(dof_handler.get_dof_handlers().at(index)),
        initialCondition<dim>(index, subset_attributes.at(index).field_type, user_inputs),
        *(solution_handler.get_solution_vector(index, dependencyType::NORMAL)));

      // TODO (landinjm): Fix so that we apply some sort of initial condition to all old
      // vector for all types.
      if (solution_handler.get_solution_vector().find(
            std::make_pair(index, dependencyType::OLD_1)) !=
          solution_handler.get_solution_vector().end())
        {
          *(solution_handler.get_solution_vector(index, dependencyType::OLD_1)) =
            *(solution_handler.get_solution_vector(index, dependencyType::NORMAL));
        }
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

PRISMS_PF_END_NAMESPACE