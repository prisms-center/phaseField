// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/user_inputs/constraint_parameters.h>
#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
BoundaryParameters::declare(dealii::ParameterHandler &parameter_handler,
                            unsigned int              n_subsections)
{
  for (unsigned int criterion_id = 0; criterion_id < n_subsections; criterion_id++)
    {
      std::string subsection_text = "boundary conditions " + std::to_string(criterion_id);
      parameter_handler.enter_subsection(subsection_text);
      {
        declare_entry(
          parameter_handler,
          "variables",
          "",
          Patterns::StringList(),
          "A list of field names that will use the constraints specified in this "
          "subsection.",
          {"fields"},
          "",
          "set variables = phi, psi");
        declare_entry(
          parameter_handler,
          "conditions",
          "Natural",
          Patterns::SelectionList("Natural|Dirichlet|Neumann|Periodic"),
          "A list of boundary conditions that correspond to the boundary ids of the "
          "mesh. For example, the first condition in the list is applied to boundary id "
          "0. If only one condition is specified, it is applied to all boundary ids.",
          {},
          "Periodic boundary conditions must be specified on the mesh! The periodic "
          "boundary conditions here are more of placeholder than anything. Additionally, "
          "if you change your mesh, your boundary ids will change.",
          "# A rectangular mesh has the following boundary id order:\n"
          "#  0 -> left (-x)\n"
          "#  1 -> right (+x)\n"
          "#  2 -> bottom (-y)\n"
          "#  3 -> top (+y)\n"
          "#  4 -> front (-z)\n"
          "#  5 -> back (+z)\n"
          "# An example set of boundary conditions in `dim = 2` could be:\n"
          "set conditions = Natural, Natural, Periodic, Periodic");
        declare_entry(parameter_handler,
                      "time dependent",
                      "false",
                      Patterns::Bool(),
                      "Whether these boundary conditions vary in time.",
                      {"time varying"},
                      "Time dependent boundary conditions will slow down simulations.");
      }
      parameter_handler.leave_subsection();
    }
}

template <unsigned int dim>
void
BoundaryParameters::assign(dealii::ParameterHandler &parameter_handler,
                           unsigned int              n_subsections)
{
  for (unsigned int criterion_id = 0; criterion_id < n_subsections; criterion_id++)
    {
      std::string subsection_text =
        "boundary conditions: " + std::to_string(criterion_id);
      parameter_handler.enter_subsection(subsection_text);
      {
        std::vector<std::string> field_names =
          dealii::Utilities::split_string_list(parameter_handler.get("variables"));
        std::vector<std::string> conditions_strings =
          dealii::Utilities::split_string_list(parameter_handler.get("conditions"));

        ComponentConditions component_conditions;
        for (unsigned int boundary_id = 0; boundary_id < conditions_strings.size();
             boundary_id++)
          {
            component_conditions[boundary_id] =
              condition_from_string(conditions_strings[boundary_id]);
          }

        // Attach conditions to fields
        for (const auto &field_comp_name : field_names)
          {
            std::string            field_name;
            std::set<unsigned int> comps;
            if (field_comp_name.ends_with(":x"))
              {
                comps      = {0};
                field_name = field_comp_name.substr(0, field_comp_name.length() - 2);
              }
            else if (field_comp_name.ends_with(":y"))
              {
                comps      = {1};
                field_name = field_comp_name.substr(0, field_comp_name.length() - 2);
              }
            else if (field_comp_name.ends_with(":z"))
              {
                comps      = {2};
                field_name = field_comp_name.substr(0, field_comp_name.length() - 2);
              }
            else
              {
                for (unsigned int comp = 0; comp < dim; ++comp)
                  {
                    comps.insert(comp);
                  }
                field_name = field_comp_name;
              }
            for (unsigned int component : comps)
              {
                if (component < dim)
                  {
                    boundary_condition_list[field_name].component_constraints[component] =
                      component_conditions;
                  }
              }
            boundary_condition_list[field_name].time_dependent =
              parameter_handler.get_bool("time dependent");
          }
      }
      parameter_handler.leave_subsection();
    }
}

void
BoundaryParameters::validate(const std::vector<FieldAttributes> &field_attributes,
                             const std::vector<SolveBlock>      &solve_blocks) const
{
  // TODO: Do this later
}

#include "user_inputs/constraint_parameters.inst"

PRISMS_PF_END_NAMESPACE
