// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/user_inputs/grain_reassignment_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
GrainReassignmentParameters::declare(dealii::ParameterHandler &parameter_handler,
                                     unsigned int              n_subsections)
{
  parameter_handler.declare_entry("use grain reassignment",
                                  "false",
                                  dealii::Patterns::Bool(),
                                  "Whether to perform grain reassignment");

  parameter_handler.declare_alias("use grain reassignment", "use_grain_reassignment");
  parameter_handler.declare_alias("use grain reassignment", "use grain_reassignment");
  parameter_handler.declare_alias("use grain reassignment", "use_grain reassignment");
  parameter_handler.declare_alias("use grain reassignment", "use grain remapping");
  parameter_handler.declare_alias("use grain reassignment", "use_grain_remapping");
  parameter_handler.declare_alias("use grain reassignment", "use_grain remapping");
  parameter_handler.declare_alias("use grain reassignment", "use grain_remapping");

  parameter_handler.enter_subsection("grain reassignment");
  {
    parameter_handler.declare_entry(
      "steps between reassignment",
      "100",
      dealii::Patterns::Integer(1),
      "The number of steps between performing grain reassignment");

    parameter_handler.declare_entry(
      "exclusion distance",
      "5.0",
      dealii::Patterns::Double(0.0),
      "The maximum distance between grains of the same order parameter");

    parameter_handler.declare_entry(
      "grain identification threshold",
      "1e-2",
      dealii::Patterns::Double(1e-15),
      "The threshold used by the flood filling algorithm for identifying grains");

    parameter_handler.declare_entry(
      "load grain structure",
      "false",
      dealii::Patterns::Bool(),
      "Whether to load the grain IDs from a file (NOT CURRENTLY IMPLEMENTED)");
  }
  parameter_handler.leave_subsection();

  // TODO: add more aliases for these parameters
}

void
GrainReassignmentParameters::assign(dealii::ParameterHandler &parameter_handler,
                                    unsigned int              n_subsections)
{
  grain_reassignment_active = parameter_handler.get_bool("use grain reassignment");

  parameter_handler.enter_subsection("grain reassignment");
  {
    reassignment_period = parameter_handler.get_integer("steps between reassignment");

    exclusion_distance = parameter_handler.get_double("exclusion distance");

    order_parameter_threshold =
      parameter_handler.get_double("grain identification threshold");

    // TODO: this currently is not implemented and does nothing
    load_grain_structure = parameter_handler.get_bool("load grain structure");
  }
  parameter_handler.leave_subsection();
}

void
GrainReassignmentParameters::validate(
  const std::vector<FieldAttributes> &field_attributes,
  const std::vector<SolveBlock>      &solve_blocks) const
{
  // TODO: Do this later
}

bool
GrainReassignmentParameters::should_perform_grain_reassignment(
  unsigned int increment) const
{
  return grain_reassignment_active && increment % reassignment_period == 0;
}

PRISMS_PF_END_NAMESPACE
