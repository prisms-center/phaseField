// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief A base class for the various subcontainers of the user inputs.
 */
template <unsigned int dim>
class ParameterBase
{
public:
  /**
   * @brief  Constructor.
   */
  ParameterBase() = default;

  /**
   * @brief Virtual destructor.
   */
  virtual ~ParameterBase() = default;

  /**
   * @brief Postprocess and validate parameters.
   */
  virtual void
  postprocess_and_validate(
    const std::map<unsigned int, VariableAttributes> &var_attributes) = 0;

  /**
   * @brief Print parameters to summary.log
   */
  virtual void
  print_parameter_summary() const
  {
    ConditionalOStreams::pout_summary()
      << "================================================\n"
      << "  Base Parameter Class - No Parameters Defined\n"
      << "================================================\n"
      << std::flush;
  };
};

PRISMS_PF_END_NAMESPACE