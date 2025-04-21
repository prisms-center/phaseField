// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/point.h>
#include <deal.II/base/types.h>
#include <deal.II/base/utilities.h>

#include <boost/algorithm/string/predicate.hpp>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

#include <cstddef>
#include <map>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Struct that stores relevant information for boundary conditions of a certain
 * field.
 */
struct boundaryCondition
{
public:
  /**
   * \brief Type of boundary condition.
   */
  enum type : std::uint8_t
  {
    UNDEFINED_BOUNDARY,
    NATURAL,
    DIRICHLET,
    PERIODIC,
    NEUMANN,
    NON_UNIFORM_DIRICHLET,
    NON_UNIFORM_NEUMANN
  };

  /**
   * \brief Test for equality of two boundary conditions.
   */
  bool
  operator==(const boundaryCondition &other) const
  {
    return boundary_condition_map == other.boundary_condition_map &&
           dirichlet_value_map == other.dirichlet_value_map;
  }

  // Map of boundary conditions and domain boundary for which they correspond to. For a
  // simple geometry like a square the boundary ids are marked, in order, by x=0, x=max,
  // y=0, y=max. More complex geometries can have somewhat arbitrary ordering, but will
  // render some of our assertions moot.
  std::map<dealii::types::boundary_id, type> boundary_condition_map;

  // A map of boundary values for dirichlet boundary conditions
  std::map<dealii::types::boundary_id, double> dirichlet_value_map;

  /**
   * \brief Enum to string for type
   */
  [[nodiscard]] std::string
  to_string(type boundary_type) const
  {
    switch (boundary_type)
      {
        case type::UNDEFINED_BOUNDARY:
          return "UNDEFINED_BOUNDARY";
        case type::NATURAL:
          return "NATURAL";
        case type::DIRICHLET:
          return "DIRICHLET";
        case type::PERIODIC:
          return "PERIODIC";
        case type::NEUMANN:
          return "NEUMANN";
        case type::NON_UNIFORM_DIRICHLET:
          return "NON_UNIFORM_DIRICHLET";
        case type::NON_UNIFORM_NEUMANN:
          return "NON_UNIFORM_NEUMANN";
        default:
          return "UNKNOWN";
      }
  }
};

/**
 * \brief Struct that holds boundary parameters.
 */
template <unsigned int dim>
struct boundaryParameters
{
public:
  using BoundaryConditionMap =
    std::map<types::index, std::map<unsigned int, boundaryCondition>>;
  using BCList         = std::map<types::index, std::map<unsigned int, std::string>>;
  using PinnedPointMap = std::map<types::index, std::pair<double, dealii::Point<dim>>>;

  /**
   * \brief Postprocess and validate parameters.
   */
  void
  postprocess_and_validate(
    const std::map<unsigned int, variableAttributes> &var_attributes);

  /**
   * \brief Check whether the boundary conditions for two fields are the same.
   */
  [[nodiscard]] bool
  check_duplicate_boundary_conditions(const types::index &index_1,
                                      const types::index &index_2) const;

  /**
   * \brief Print parameters to summary.log
   */
  void
  print_parameter_summary() const;

  // Map of unfiltered boundary conditions strings. The first key is the global index. The
  // second key is the number of dimensions.
  BCList BC_list;

  // Map of pinned points. The first key is the global index. The pair is the pinned
  // value and point.
  PinnedPointMap pinned_point_list = {};

  // Map of boundary conditions. The first key is the global index. The second key is the
  // number of dimensions.
  BoundaryConditionMap boundary_condition_list;

private:
  /**
   * \brief Set the boundary for a single component of a field index.
   */
  void
  set_boundary(const std::string  &BC_string,
               const types::index &index,
               const unsigned int &component);

  /**
   * \brief Perform a check on the boundary conditions to ensure that they are valid
   */
  void
  validate_boundary_conditions() const;
};

template <unsigned int dim>
inline void
boundaryParameters<dim>::postprocess_and_validate(
  const std::map<unsigned int, variableAttributes> &var_attributes)
{
  for (const auto &[index, variable] : var_attributes)
    {
      // Ensure that boundary conditions are specified for all variables and their
      // components
      if (variable.field_type == fieldType::VECTOR)
        {
          for (unsigned int i = 0; i < dim; i++)
            {
              if (variable.is_postprocess)
                {
                  set_boundary("NATURAL", variable.field_index, i);
                }
              else
                {
                  AssertThrow(BC_list.contains(variable.field_index),
                              dealii::ExcMessage("Invalid entry"));
                  AssertThrow(BC_list.at(variable.field_index).contains(i),
                              dealii::ExcMessage("Invalid entry"));
                  AssertThrow(!BC_list.at(variable.field_index).at(i).empty(),
                              dealii::ExcMessage(
                                "Boundary conditions must be specified "
                                "for all components in all vector field."));

                  set_boundary(BC_list.at(variable.field_index).at(i),
                               variable.field_index,
                               i);
                }
            }
        }
      else
        {
          if (variable.is_postprocess)
            {
              set_boundary("NATURAL", variable.field_index, 0);
            }
          else
            {
              AssertThrow(BC_list.contains(variable.field_index),
                          dealii::ExcMessage("Invalid entry"));
              AssertThrow(BC_list.at(variable.field_index).contains(0),
                          dealii::ExcMessage("Invalid entry"));
              AssertThrow(!BC_list.at(variable.field_index).at(0).empty(),
                          dealii::ExcMessage("Boundary conditions must be specified "
                                             "for all scalar fields."));

              set_boundary(BC_list.at(variable.field_index).at(0),
                           variable.field_index,
                           0);
            }
        }
    }

  // Validate boundary conditions
  validate_boundary_conditions();

  // Clear the BC_list now that it's no longer necessary
  BC_list.clear();

#ifdef ADDITIONAL_OPTIMIZATIONS
  // Check if any fields are duplicates in terms of boundary conditions
  // TODO (landinjm): Clean this up
  for (const auto &[index_1, variable_1] : var_attributes)
    {
      // Skip is the duplicate index has already been assigned
      if (variable_1.duplicate_field_index != numbers::invalid_index)
        {
          continue;
        }
      if (variable_1.is_postprocess)
        {
          continue;
        }

      const auto field_type_1 = variable_1.field_type;

      for (const auto &[index_2, variable_2] : var_attributes)
        {
          if (variable_2.is_postprocess)
            {
              continue;
            }

          bool is_duplicate = false;

          const auto field_type_2 = variable_2.field_type;

          is_duplicate = field_type_1 == field_type_2 &&
                         check_duplicate_boundary_conditions(index_1, index_2);

          if (is_duplicate)
            {
              conditionalOStreams::pout_verbose()
                << "Field " << variable_1.name << " has the same boundary conditions as "
                << variable_2.name << ". Using optimizations...\n";
              variable_2.duplicate_field_index = index_1;
            }
        }
    }
#endif
}

template <unsigned int dim>
inline bool
boundaryParameters<dim>::check_duplicate_boundary_conditions(
  const types::index &index_1,
  const types::index &index_2) const
{
  // If the indices are the same return false
  if (index_1 == index_2)
    {
      return false;
    }

  Assert(boundary_condition_list.contains(index_1),
         dealii::ExcMessage("Invalid entry for index = " + std::to_string(index_1)));
  Assert(boundary_condition_list.contains(index_2),
         dealii::ExcMessage("Invalid entry for index = " + std::to_string(index_2)));

  bool is_duplicate = false;

  // First check the boundary_condition_list
  const auto &boundary_condition_1 = boundary_condition_list.at(index_1);
  const auto &boundary_condition_2 = boundary_condition_list.at(index_2);

  is_duplicate = boundary_condition_1 == boundary_condition_2;

  // Check the pinned points
  if (pinned_point_list.contains(index_1) && pinned_point_list.contains(index_2))
    {
      is_duplicate =
        is_duplicate && pinned_point_list.at(index_1) == pinned_point_list.at(index_2);
    }

  return is_duplicate;
}

template <unsigned int dim>
inline void
boundaryParameters<dim>::print_parameter_summary() const
{
  conditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Boundary Parameters\n"
    << "================================================\n";

  for (const auto &[index, component_map] : boundary_condition_list)
    {
      conditionalOStreams::pout_summary() << "Index: " << index << "\n";
      for (const auto &[component, boundary_condition] : component_map)
        {
          conditionalOStreams::pout_summary() << "  Component: " << component << "\n";
          for (const auto &[domain_id, boundary_type] :
               boundary_condition.boundary_condition_map)
            {
              conditionalOStreams::pout_summary()
                << "    Boundary id: " << domain_id << "    "
                << "Condition: " << boundary_condition.to_string(boundary_type);
              if (boundary_type == boundaryCondition::type::DIRICHLET)
                {
                  conditionalOStreams::pout_summary()
                    << " = " << boundary_condition.dirichlet_value_map.at(domain_id);
                }
              conditionalOStreams::pout_summary() << "\n";
            }
        }
    }

  if (!pinned_point_list.empty())
    {
      conditionalOStreams::pout_summary() << "Pinned field index: ";
    }
  for (const auto &[index, point_value_map] : pinned_point_list)
    {
      conditionalOStreams::pout_summary()
        << index << "\n"
        << "  Value: " << point_value_map.first << "\n"
        << "  Point: " << point_value_map.second << "\n";
    }
  conditionalOStreams::pout_summary() << "\n" << std::flush;
}

template <unsigned int dim>
inline void
boundaryParameters<dim>::set_boundary(const std::string  &BC_string,
                                      const types::index &index,
                                      const unsigned int &component)
{
  // Split string
  auto BC_string_list = dealii::Utilities::split_string_list(BC_string);

  // Check that there is either 1 or 2*dim entries in the vector. This can be changed
  // later to support other geometries.
  AssertThrow(BC_string_list.size() == 1 ||
                BC_string_list.size() == static_cast<std::size_t>(2 * dim),
              dealii::ExcMessage(
                "Either 1 or 2*dim boundary conditions must be specified."));

  // If there is only 1 boundary condition resize BC_string_list, copying the first
  // entry.
  if (BC_string_list.size() == 1)
    {
      BC_string_list.resize(static_cast<std::size_t>(2 * dim), BC_string_list[0]);
    }

  // Assign boundary condition
  boundaryCondition condition;
  for (unsigned int i = 0; i < (2 * dim); i++)
    {
      const std::string dirichlet = "DIRICHLET";
      const std::string neumann   = "NEUMANN";

      if (boost::iequals(BC_string_list[i], "NATURAL"))
        {
          condition.boundary_condition_map.emplace(i, boundaryCondition::type::NATURAL);
        }
      else if (boost::iequals(BC_string_list[i].substr(0, dirichlet.size()), dirichlet))
        {
          condition.boundary_condition_map.emplace(i, boundaryCondition::type::DIRICHLET);
          std::string dirichlet_value =
            BC_string_list[i].substr(dirichlet.size() + 1, BC_string_list[i].size());
          dirichlet_value = dealii::Utilities::trim(dirichlet_value);
          condition.dirichlet_value_map.emplace(i,
                                                dealii::Utilities::string_to_double(
                                                  dirichlet_value));
        }
      else if (boost::iequals(BC_string_list[i], "PERIODIC"))
        {
          condition.boundary_condition_map.emplace(i, boundaryCondition::type::PERIODIC);
        }
      else if (boost::iequals(BC_string_list[i].substr(0, neumann.size()), neumann))
        {
          AssertThrow(false, FeatureNotImplemented("Neumann boundary conditions"));
        }
      else if (boost::iequals(BC_string_list[i], "NON_UNIFORM_DIRICHLET"))
        {
          condition.boundary_condition_map
            .emplace(i, boundaryCondition::type::NON_UNIFORM_DIRICHLET);
        }
      else if (boost::iequals(BC_string_list[i], "NON_UNIFORM_NEUMANN"))
        {
          AssertThrow(false,
                      FeatureNotImplemented("Nonuniform neumann boundary conditions"));
        }
      else
        {
          AssertThrow(false,
                      dealii::ExcMessage("Invalid boundary condition " +
                                         BC_string_list[i]));
        }
      // If periodic boundary conditions are used, ensure that they are applied on
      // both sides of the domain.
      if (i % 2 == 0)
        {
          AssertThrow(boost::iequals(BC_string_list[i], "PERIODIC") ==
                        boost::iequals(BC_string_list[i + 1], "PERIODIC"),
                      dealii::ExcMessage("Periodic boundary condition must be "
                                         "specified on both sides of domain"));
        }
    }

  boundary_condition_list[index].emplace(component, condition);
}

template <unsigned int dim>
inline void
boundaryParameters<dim>::validate_boundary_conditions() const
{
  // Throw a warning if the pinned point is not on a vertex
  // TODO (landinjm): This should be fixed
  for (const auto &[index, point_value_map] : pinned_point_list)
    {
      const auto               point = point_value_map.second;
      const dealii::Point<dim> origin {};
      AssertThrow(point == origin,
                  dealii::ExcMessage("Pinned point must be on the origin"));
    }

  // Throw a warning if only some fields have periodic boundary conditions
  std::vector<bool> periodic_ids(static_cast<dealii::types::boundary_id>(2 * dim), false);
  for (const auto &[index, component_map] : boundary_condition_list)
    {
      for (const auto &[component, boundary_condition] : component_map)
        {
          for (const auto &[domain_id, boundary_type] :
               boundary_condition.boundary_condition_map)
            {
              if (boundary_type == boundaryCondition::type::PERIODIC)
                {
                  periodic_ids[domain_id] = true;
                }
            }
        }
    }
  for (const auto &[index, component_map] : boundary_condition_list)
    {
      for (const auto &[component, boundary_condition] : component_map)
        {
          for (const auto &[domain_id, boundary_type] :
               boundary_condition.boundary_condition_map)
            {
              if (boundary_type != boundaryCondition::type::PERIODIC &&
                  periodic_ids[domain_id])
                {
                  AssertThrow(
                    false,
                    dealii::ExcMessage(
                      "All fields for a given domain id (side) must have periodic "
                      "boundary conditions if any field has periodic boundary "
                      "conditions"));
                }
            }
        }
    }
}

PRISMS_PF_END_NAMESPACE