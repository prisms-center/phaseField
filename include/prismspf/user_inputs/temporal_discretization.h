// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef temporal_discretization_h
#define temporal_discretization_h

#include <prismspf/config.h>
#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/variable_attributes.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Struct that holds temporal discretization parameters.
 */
struct temporalDiscretization
{
public:
  /**
   * \brief Postprocess and validate parameters.
   */
  void
  postprocess_and_validate(
    const std::map<unsigned int, variableAttributes> &var_attributes);

  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  // Final time
  double final_time = 0.0;

  // Total number of increments
  unsigned int total_increments = 0;

  // Timestep
  mutable double dt = 0.0;

  // The current increment
  mutable unsigned int increment = 0;

  // The current time
  mutable double time = 0.0;
};

inline void
temporalDiscretization::postprocess_and_validate(
  const std::map<unsigned int, variableAttributes> &var_attributes)
{
  // If all of the variables are `TIME_INDEPENDENT`, `AUXILIARY`, or `CONSTANT` then
  // total_increments should be 1 and final_time should be 0
  bool only_time_independent_pdes = true;
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.pde_type == PDEType::EXPLICIT_TIME_DEPENDENT ||
          variable.pde_type == PDEType::IMPLICIT_TIME_DEPENDENT)
        {
          only_time_independent_pdes = false;
          break;
        }
    }
  if (only_time_independent_pdes)
    {
      total_increments = 1;
      return;
    }

  // Check that the timestep is greater than zero
  AssertThrow(dt > 0.0,
              dealii::ExcMessage(
                "The timestep must be greater than zero for transient problems."));

  // Pick the maximum specified time since the default values are zero
  final_time       = std::max(final_time, dt * total_increments);
  total_increments = static_cast<unsigned int>(std::ceil(final_time / dt));
}

inline void
temporalDiscretization::print_parameter_summary() const
{
  conditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Temporal Discretization\n"
    << "================================================\n"
    << "Timestep: " << dt << "\n"
    << "Total increments: " << total_increments << "\n"
    << "Final time: " << final_time << "\n\n"
    << std::flush;
}

PRISMS_PF_END_NAMESPACE

#endif