// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <boost/algorithm/string/predicate.hpp>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/config.h>

#include <algorithm>
#include <map>
#include <string>
#include <unordered_set>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Condition of boundary condition.
 */
enum Condition : std::uint8_t
{
  Natural,
  Dirichlet,
  Neumann,
  TimeDependentDirichlet,
  TimeDependentNeumann,
  UniformDirichlet,
  UniformNeumann,
  Periodic,
};

/**
 * @brief Enum to string for type
 */
[[nodiscard]] inline std::string
to_string(Condition boundary_type)
{
  switch (boundary_type)
    {
      case Condition::Natural:
        return "Natural";
      case Condition::Dirichlet:
        return "Dirichlet";
      case Condition::Neumann:
        return "Neumann";
      case Condition::TimeDependentDirichlet:
        return "TimeDependentDirichlet";
      case Condition::TimeDependentNeumann:
        return "TimeDependentNeumann";
      case Condition::UniformDirichlet:
        return "UniformDirichlet";
      case Condition::UniformNeumann:
        return "UniformNeumann";
      case Condition::Periodic:
        return "Periodic";
      default:
        return "UNKNOWN";
    }
}

/**
 * @brief Enum to string for type
 */
[[nodiscard]] inline Condition
condition_from_string(const std::string &boundary_string)
{
  if (boundary_string == "Natural")
    {
      return Condition::Natural;
    }
  if (boundary_string == "Dirichlet")
    {
      return Condition::Dirichlet;
    }
  if (boundary_string == "Neumann")
    {
      return Condition::Neumann;
    }
  if (boundary_string == "TimeDependentDirichlet")
    {
      return Condition::TimeDependentDirichlet;
    }
  if (boundary_string == "TimeDependentNeumann")
    {
      return Condition::TimeDependentNeumann;
    }
  if (boundary_string == "UniformDirichlet")
    {
      return Condition::UniformDirichlet;
    }
  if (boundary_string == "UniformNeumann")
    {
      return Condition::UniformNeumann;
    }
  if (boundary_string == "Periodic")
    {
      return Condition::Periodic;
    }
  AssertThrow(false, dealii::ExcMessage("Invalid boundary condition " + boundary_string));
  return Condition::Natural;
}

/**
 * @brief Struct that stores relevant information for boundary conditions of a certain
 * field.
 */
struct ComponentConditions
{
  // Map of boundary conditions and domain boundary for which they correspond to. For a
  // simple geometry like a square the boundary ids are marked, in order, by x=0, x=max,
  // y=0, y=max. More complex geometries can have somewhat arbitrary ordering, but will
  // render some of our assertions moot.
  std::map<unsigned int, Condition> conditions;

  [[nodiscard]] bool
  has_time_dependent_bcs() const
  {
    return std::any_of(conditions.begin(),
                       conditions.end(),
                       [](const auto &dir_cond)
                       {
                         return dir_cond.second == Condition::TimeDependentDirichlet ||
                                dir_cond.second == Condition::TimeDependentNeumann;
                       });
  }
};

template <unsigned int dim>
struct FieldConstraints
{
  std::array<ComponentConditions, dim> component_constraints;

  // std::unordered_set<std::pair<dealii::Point<dim>, dealii::Tensor<1, dim>>>
  // pinned_points;

  [[nodiscard]] bool
  has_time_dependent_bcs() const
  {
    return std::any_of(component_constraints.begin(),
                       component_constraints.end(),
                       [](const ComponentConditions &comp)
                       {
                         return comp.has_time_dependent_bcs();
                       });
  }
};

/**
 * @brief Struct that holds boundary parameters.
 */
template <unsigned int dim>
struct BoundaryParameters
{
public:
  /**
   * @brief Postprocess and validate parameters.
   */
  void
  validate();

  /**
   * @brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  /**
   * @brief Whether there are time-dependent boundary conditions.
   */
  [[nodiscard]] bool
  has_time_dependent_bcs() const
  { // todo
    return false;
  }

  /**
   * @brief Declare the parameters to be read from an input file.
   */
  void
  declare_parameters(dealii::ParameterHandler &parameter_handler,
                     unsigned int              max_criteria = 5) const;

  /**
   * @brief Assign the parameters read from an input file to this object.
   */
  void
  assign_parameters(dealii::ParameterHandler &parameter_handler,
                    unsigned int              max_criteria = 5);

  // Map of boundary conditions. The first key is the field index.
  std::unordered_map<std::string, FieldConstraints<dim>> boundary_condition_list;
};

template <unsigned int dim>
inline void
BoundaryParameters<dim>::validate()
{ // todo
}

template <unsigned int dim>
inline void
BoundaryParameters<dim>::print_parameter_summary() const
{
  ConditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Boundary Parameters\n"
    << "================================================\n";

  for (const auto &[index, component_map] : boundary_condition_list)
    {
      // Todo
    }
  // TODO print pinned points

  ConditionalOStreams::pout_summary() << "\n" << std::flush;
}

template <unsigned int dim>
inline void
BoundaryParameters<dim>::declare_parameters(dealii::ParameterHandler &parameter_handler,
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
          dealii::Patterns::Anything(),
          "The names of the fields that will use these constraints.");
        parameter_handler.declare_entry("conditions",
                                        "",
                                        dealii::Patterns::Anything(),
                                        "List of conditions.");
        parameter_handler.enter_subsection("pinning point");
        {
          parameter_handler.declare_entry("enable pinned point",
                                          "false",
                                          dealii::Patterns::Bool(),
                                          "Whether to use a pinned point");

          parameter_handler.declare_entry("x",
                                          "0.0",
                                          dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                          "X-coordinate of the point");
          parameter_handler.declare_entry("y",
                                          "0.0",
                                          dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                          "Y-coordinate of the point");
          parameter_handler.declare_entry("z",
                                          "0.0",
                                          dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                          "Z-coordinate of the point");

          parameter_handler.declare_entry("value",
                                          "2147483647",
                                          dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                          "Value of pinned point.");

          parameter_handler.declare_entry("x value",
                                          "2147483647",
                                          dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                          "Value of pinned point for the x-component.");
          parameter_handler.declare_entry("y value",
                                          "2147483647",
                                          dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                          "Value of pinned point for the y-component.");
          parameter_handler.declare_entry("z value",
                                          "2147483647",
                                          dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                          "Value of pinned point for the z-component.");
        }
        parameter_handler.leave_subsection();
      }
      parameter_handler.leave_subsection();
    }
}

template <unsigned int dim>
inline void
BoundaryParameters<dim>::assign_parameters(dealii::ParameterHandler &parameter_handler,
                                           unsigned int              max_criteria)
{
  static const std::vector<std::string> axis_labels = {"x", "y", "z"};
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
        if (conditions_strings.size() == 1)
          {
            // all the same
            for (unsigned int boundary_id = 0; boundary_id < 2 * dim; boundary_id++)
              {
                component_conditions.conditions[boundary_id] =
                  condition_from_string(conditions_strings[0]);
              }
          }
        else
          {
            for (unsigned int boundary_id = 0; boundary_id < conditions_strings.size();
                 boundary_id++)
              {
                component_conditions.conditions[boundary_id] =
                  condition_from_string(conditions_strings[boundary_id]);
              }
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

PRISMS_PF_END_NAMESPACE
