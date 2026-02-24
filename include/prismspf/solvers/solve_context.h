// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/constraint_manager.h>
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/pde_operator_base.h>
#include <prismspf/core/simulation_timer.h>
#include <prismspf/core/solution_indexer.h>
#include <prismspf/core/triangulation_manager.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include "prismspf/core/field_attributes.h"
#include "prismspf/core/invm_manager.h"

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class provides context for a solver with ptrs to all the relevant
 * dependencies.
 *
 * The context in this case refers to all the finite element machinery needed to solve the
 * fields. For example, it contains the triangulation manager and pde operator that
 * evaluates the user-specified PDEs.
 */
template <unsigned int dim, unsigned int degree, typename number>
class SolveContext
{
public:
  /**
   * @brief Constructor.
   */
  SolveContext(std::vector<FieldAttributes>            _field_attributes,
               const UserInputParameters<dim>         &_user_inputs,
               TriangulationManager<dim>              &_triangulation_manager,
               DofManager<dim>                        &_dof_manager,
               ConstraintManager<dim, degree, number> &_constraint_manager,
               SolutionIndexer<dim, number>           &_solution_indexer,
               std::shared_ptr<PDEOperatorBase<dim, degree, number>> _pde_operator)
    : field_attributes(std::move(_field_attributes))
    , user_inputs(&_user_inputs)
    , triangulation_manager(&_triangulation_manager)
    , dof_manager(&_dof_manager)
    , constraint_manager(&_constraint_manager)
    , solution_indexer(&_solution_indexer)
    , invm_manager(*dof_manager, true, true)
    , sim_timer(user_inputs->get_temporal_discretization().dt)
    , pde_operator(_pde_operator) {};

  /**
   * @brief Get the field attributes.
   */
  [[nodiscard]] const std::vector<FieldAttributes> &
  get_field_attributes() const
  {
    return field_attributes;
  }

  /**
   * @brief Get the user-inputs.
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const
  {
    Assert(user_inputs != nullptr, dealii::ExcNotInitialized());
    return *user_inputs;
  }

  /**
   * @brief Get the triangulation manager.
   */
  [[nodiscard]] const TriangulationManager<dim> &
  get_triangulation_manager() const
  {
    Assert(triangulation_manager != nullptr, dealii::ExcNotInitialized());
    return *triangulation_manager;
  }

  /**
   * @brief Get the triangulation manager.
   */
  [[nodiscard]] TriangulationManager<dim> &
  get_triangulation_manager()
  {
    Assert(triangulation_manager != nullptr, dealii::ExcNotInitialized());
    return *triangulation_manager;
  }

  /**
   * @brief Get the dof manager.
   */
  [[nodiscard]] const DofManager<dim> &
  get_dof_manager() const
  {
    Assert(dof_manager != nullptr, dealii::ExcNotInitialized());
    return *dof_manager;
  }

  /**
   * @brief Get the dof manager.
   */
  [[nodiscard]] DofManager<dim> &
  get_dof_manager()
  {
    Assert(dof_manager != nullptr, dealii::ExcNotInitialized());
    return *dof_manager;
  }

  /**
   * @brief Get the constraint manager.
   */
  [[nodiscard]] const ConstraintManager<dim, degree, number> &
  get_constraint_manager() const
  {
    Assert(constraint_manager != nullptr, dealii::ExcNotInitialized());
    return *constraint_manager;
  }

  /**
   * @brief Get the constraint manager.
   */
  [[nodiscard]] ConstraintManager<dim, degree, number> &
  get_constraint_manager()
  {
    Assert(constraint_manager != nullptr, dealii::ExcNotInitialized());
    return *constraint_manager;
  }

  /**
   * @brief Get the solution manager.
   */
  [[nodiscard]] SolutionIndexer<dim, number> &
  get_solution_indexer() const
  {
    Assert(solution_indexer != nullptr, dealii::ExcNotInitialized());
    return *solution_indexer;
  }

  /**
   * @brief Get the invm manager.
   */
  [[nodiscard]] const InvMManager<dim, degree, number> &
  get_invm_manager() const
  {
    return invm_manager;
  }

  /**
   * @brief Get the invm manager.
   */
  [[nodiscard]] InvMManager<dim, degree, number> &
  get_invm_manager()
  {
    return invm_manager;
  }

  /**
   * @brief Get the simulation timer.
   */
  [[nodiscard]] const SimulationTimer &
  get_simulation_timer() const
  {
    return sim_timer;
  }

  /**
   * @brief Get the simulation timer.
   */
  [[nodiscard]] SimulationTimer &
  get_simulation_timer()
  {
    return sim_timer;
  }

  /**
   * @brief Get a shared pointer to the pde operator.
   */
  [[nodiscard]] const std::shared_ptr<PDEOperatorBase<dim, degree, number>> &
  get_pde_operator() const
  {
    Assert(pde_operator != nullptr, dealii::ExcNotInitialized());
    return pde_operator;
  }

private:
  /**
   * @brief Field attributes.
   */
  std::vector<FieldAttributes> field_attributes;

  /**
   * @brief User-inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * @brief Triangulation manager.
   */
  TriangulationManager<dim> *triangulation_manager;

  /**
   * @brief DoF manager.
   */
  DofManager<dim> *dof_manager;

  /**
   * @brief Constraint manager.
   */
  ConstraintManager<dim, degree, number> *constraint_manager;

  /**
   * @brief Solution manager.
   */
  SolutionIndexer<dim, number> *solution_indexer;

  /**
   * @brief Solution manager.
   */
  InvMManager<dim, degree, number> invm_manager;

  /**
   * @brief Simulation timer.
   */
  SimulationTimer sim_timer;

  /**
   * @brief PDE operator.
   */
  std::shared_ptr<PDEOperatorBase<dim, degree, number>> pde_operator;
};

PRISMS_PF_END_NAMESPACE
