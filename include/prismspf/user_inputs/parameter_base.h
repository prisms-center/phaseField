// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <ostream>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief A base class for the various subcontainers of the user inputs.
 *
 * There are five important steps that all child classes must achieve.
 * 1. Declare entries in the `dealii::ParameterHandler`. We can preprocess the data so
 * users don't give us nonsensical values.
 * 2. Read entries from the `dealii:ParameterHandler` into the class.
 * 3. Postprocess the data into something more useful for PRISMS-PF since we know what
 * fields exist.
 * 4. Validate the parameters to make sure everything is correct.
 * 5. Package the parameters into a minimal struct so this class can be deconstructed.
 */
class ParameterBase
{
public:
  /**
   * @brief Constructor.
   */
  ParameterBase() = default;

  /**
   * @brief Virtual destructor.
   */
  virtual ~ParameterBase() = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so we don't copy this base.
   */
  ParameterBase(const ParameterBase &) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so we don't copy this base.
   */
  ParameterBase &
  operator=(const ParameterBase &) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so we don't move the base.
   */
  ParameterBase(ParameterBase &&) = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so we don't move the base.
   */
  ParameterBase &
  operator=(ParameterBase &&) = delete;

  /**
   * @brief Declare parameters.
   */
  virtual void
  declare(dealii::ParameterHandler &parameter_handler) const = 0;

  /**
   * @brief Read parameters.
   */
  virtual void
  read(dealii::ParameterHandler &parameter_handler) = 0;

  /**
   * @brief Validate parameters.
   */
  virtual void
  validate(const std::vector<VariableAttributes> &field_attributes) const = 0;

  /**
   * @brief Print parameters to summary.log
   */
  virtual void
  print_parameter_summary() const
  {
    ConditionalOStreams::pout_summary()
      << "================================================\n"
      << "  Base Parameter Class - Parameters Undefined\n"
      << "================================================\n"
      << std::flush;
  };
};

PRISMS_PF_END_NAMESPACE
