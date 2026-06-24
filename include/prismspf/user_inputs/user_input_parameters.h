// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solve_block.h>

#include <prismspf/user_inputs/checkpoint_parameters.h>
#include <prismspf/user_inputs/constraint_parameters.h>
#include <prismspf/user_inputs/linear_solve_parameters.h>
#include <prismspf/user_inputs/load_initial_condition_parameters.h>
#include <prismspf/user_inputs/miscellaneous_parameters.h>
#include <prismspf/user_inputs/nonlinear_solve_parameters.h>
#include <prismspf/user_inputs/nucleation_parameters.h>
#include <prismspf/user_inputs/output_parameters.h>
#include <prismspf/user_inputs/spatial_discretization.h>
#include <prismspf/user_inputs/temporal_discretization.h>
#include <prismspf/user_inputs/user_constants.h>

#include <prismspf/config.h>

#include "prismspf/core/types.h"

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * Helper function for declaring multiple aliases
 */
void
declare_aliases(dealii::ParameterHandler       &parameter_handler,
                const std::string              &existing_entry_name,
                const std::vector<std::string> &aliases)
{
  for (const auto &alias : aliases)
    {
      parameter_handler.declare_alias(existing_entry_name, alias);
    }
}

/**
 * @brief Virtual base class for parameter groups.
 *
 * This virtual base class ensures that we always have the following execution order.
 * 1. predeclare
 * 2. preassign
 * 3. declare
 * 4. assign
 * 5. validate
 */
struct ParameterBase
{
  virtual ~ParameterBase() = default;

  /**
   * @brief Declare the parameters to be read from file.
   *
   * Unlike `declare` this step comes first so that we can read certain parameters and
   * things easier later on. For example, there's the mesh type parameter. We must first
   * read the value of this parameter before we generate the corresponding mesh object and
   * its parameters.
   */
  virtual void
  predeclare(dealii::ParameterHandler &parameter_handler) const = 0;

  /**
   * @brief Assign the parameters from file.
   */
  virtual void
  preassign(dealii::ParameterHandler &parameter_handler) = 0;

  /**
   * @brief Declare the parameters to be read from file.
   */
  virtual void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int              max_criteria = Numbers::max_subsections) const = 0;

  /**
   * @brief Assign the parameters from file.
   */
  virtual void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections) = 0;

  /**
   * @brief Validate.
   */
  virtual void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const = 0;
};

template <unsigned int dim>
struct UserInputParameters : public ParameterBase
{
public:
  /**
   * @brief Constructor.
   */
  UserInputParameters() = default;

  /**
   * @brief Constructor.
   *
   * Read in user input parameters from file and load them into member variables.
   */
  explicit UserInputParameters(const std::string &file_name,
                               unsigned int max_criteria = Numbers::max_subsections) {};

  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  predeclare(dealii::ParameterHandler &parameter_handler) const override {};

  /**
   * @brief Assign the parameters from file.
   */
  void
  preassign(dealii::ParameterHandler &parameter_handler) override {};

  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override {};

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections) override {};

  /**
   * @brief Ensure that the parameters are compatible with a set of fields and solvers.
   *
   * This is an optional step.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override
  {
    spatial_discretization.validate(field_attributes, solve_blocks);
  }

  SpatialDiscretization<dim> spatial_discretization;
  TemporalDiscretization     temporal_discretization;

  LinearSolveParameters          linear_solve_parameters;
  NonlinearSolveParameterSet     nonlinear_solve_parameters;
  OutputParameters               output_parameters;
  CheckpointParameters           checkpoint_parameters;
  BoundaryParameters<dim>        boundary_parameters;
  LoadInitialConditionParameters load_ic_parameters;
  NucleationParameters           nucleation_parameters;
  MiscellaneousParameters        misc_parameters;
  UserConstants<dim>             user_constants;
};

template <unsigned int dim>
inline UserInputParameters<dim>::UserInputParameters(const std::string &file_name,
                                                     unsigned int       max_criteria)
{
  // user_constants is a special little princess that needs to know the file contents
  // before we declare the parameters.
  user_constants.file_name = file_name;
  dealii::ParameterHandler parameter_handler;
  declare_parameters(parameter_handler, max_criteria);
  parameter_handler.parse_input(file_name);
  assign_parameters(parameter_handler, max_criteria);
}

PRISMS_PF_END_NAMESPACE
