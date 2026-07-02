#include <prismspf/user_inputs/constraint_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

bool
ComponentConditions::has_time_dependent_bcs() const
{
  return std::any_of(conditions.begin(),
                     conditions.end(),
                     [](const auto &_condition)
                     {
                       return _condition.second == Condition::TimeDependentDirichlet ||
                              _condition.second == Condition::TimeDependentNeumann;
                     });
}

template <unsigned int dim>
bool
FieldConstraints<dim>::has_time_dependent_bcs() const
{
  return std::any_of(component_constraints.begin(),
                     component_constraints.end(),
                     [](const ComponentConditions &component)
                     {
                       return component.has_time_dependent_bcs();
                     });
}

template <unsigned int dim>
void
BoundaryParameters<dim>::predeclare(dealii::ParameterHandler &parameter_handler) const
{
  AssertThrow(false, dealii::ExcNotImplemented());
}

template <unsigned int dim>
void
BoundaryParameters<dim>::preassign(dealii::ParameterHandler &parameter_handler)
{
  AssertThrow(false, dealii::ExcNotImplemented());
}

template <unsigned int dim>
void
BoundaryParameters<dim>::declare(dealii::ParameterHandler &parameter_handler,
                                 unsigned int              max_criteria) const
{
  for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
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
      }
      parameter_handler.leave_subsection();
    }
}

template <unsigned int dim>
void
BoundaryParameters<dim>::assign(dealii::ParameterHandler &parameter_handler,
                                unsigned int              max_criteria)
{
  for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
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
            component_conditions.conditions[boundary_id] =
              condition_from_string(conditions_strings[boundary_id]);
          }

        // Attach conditions to fields
        for (const auto &field_comp_name : field_names)
          {
            int                    pos = field_comp_name.length() - 2;
            const std::string      end = field_comp_name.substr(pos > 0 ? pos : 0);
            std::string            field_name;
            std::set<unsigned int> comps;
            if (end == ":x")
              {
                comps      = {0};
                field_name = field_comp_name.substr(0, pos);
              }
            else if (end == ":y")
              {
                comps      = {1};
                field_name = field_comp_name.substr(0, pos);
              }
            else if (end == ":z")
              {
                comps      = {2};
                field_name = field_comp_name.substr(0, pos);
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
                    boundary_condition_list[field_name].component_constraints.at(
                      component) = component_conditions;
                  }
              }
          }
      }
      parameter_handler.leave_subsection();
    }
}

template <unsigned int dim>
void
BoundaryParameters<dim>::validate(const std::vector<FieldAttributes> &field_attributes,
                                  const std::vector<SolveBlock>      &solve_blocks) const
{
  // TODO: Do this later
}

PRISMS_PF_END_NAMESPACE
