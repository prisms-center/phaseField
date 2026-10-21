// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/user_inputs/grain_reassignment_parameters.h>
#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
GrainReassignmentParameters::declare(dealii::ParameterHandler &parameter_handler,
                                     unsigned int              n_subsections)
{
  declare_entry(parameter_handler,
                "use grain reassignment",
                "false",
                Patterns::Bool(),
                "Whether to perform grain reassignment.",
                {"use grain remapping"});

  parameter_handler.enter_subsection("grain reassignment");
  {
    declare_entry(parameter_handler,
                  "reassignment stride",
                  "2147483647",
                  Patterns::PositiveInteger(),
                  "The number of steps between grain reassignment attempts.",
                  {"steps between reassignment",
                   "number of steps between reassignment",
                   "num steps between reassignment",
                   "reassignment period"});
    declare_entry(parameter_handler,
                  "exclusion distance",
                  "0.0",
                  Patterns::PositiveDouble(),
                  "The minimum distance between grains of the same order parameter.");
    declare_entry(
      parameter_handler,
      "grain identification threshold",
      "0.0",
      Patterns::PositiveDouble(),
      "The threshold used by the flood filling algorithm for identifying grains.",
      {"identification threshold"});
    declare_entry(parameter_handler,
                  "load grain structure",
                  "false",
                  Patterns::Bool(),
                  "Whether to load the grain IDs from a file.",
                  {"load structure"});
  }
  parameter_handler.leave_subsection();
}

void
GrainReassignmentParameters::assign(dealii::ParameterHandler &parameter_handler,
                                    unsigned int              n_subsections)
{
  grain_reassignment_active = parameter_handler.get_bool("use grain reassignment");

  parameter_handler.enter_subsection("grain reassignment");
  {
    reassignment_period = parameter_handler.get_integer("reassignment stride");

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
