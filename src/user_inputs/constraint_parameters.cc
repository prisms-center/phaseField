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
      std::string subsection_text =
        "boundary conditions: " + std::to_string(criterion_id);
      parameter_handler.enter_subsection(subsection_text);
      {
        parameter_handler.declare_entry(
          "variables",
          "",
          dealii::Patterns::List(dealii::Patterns::Anything(), 0, INT_MAX, ","),
          "The names of the fields that will use these constraints.");
        parameter_handler.declare_entry(
          "conditions",
          "",
          dealii::Patterns::List(dealii::Patterns::Anything(), 0, INT_MAX, ","),
          "List of conditions.");
        parameter_handler.declare_entry("time dependent",
                                        "false",
                                        dealii::Patterns::Bool(),
                                        "Whether these conditions vary in time.");
        declare_aliases(parameter_handler, "time dependent");
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
