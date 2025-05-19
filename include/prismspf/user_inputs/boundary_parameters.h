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
#include <variant>
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
    NON_UNIFORM_NEUMANN,
    TIME_DEPENDENT_NON_UNIFORM_DIRICHLET,
    TIME_DEPENDENT_NON_UNIFORM_NEUMANN,
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

  /**
   * \brief Get the map of boundary conditions.
   */
  [[nodiscard]] const std::map<dealii::types::boundary_id, type> &
  get_boundary_condition_map() const
  {
    return boundary_condition_map;
  }

  /**
   * \brief Add a boundary conditions.
   */
  void
  add_boundary_condition(dealii::types::boundary_id boundary_id, type boundary_type)
  {
    boundary_condition_map.emplace(boundary_id, boundary_type);
  }

  /**
   * \brief Get the value for a homogenous dirichlet boundary condition.
   */
  [[nodiscard]] double
  get_dirichlet_value(dealii::types::boundary_id boundary_id) const
  {
    Assert(boundary_condition_map.at(boundary_id) == type::DIRICHLET,
           dealii::ExcMessage("Boundary condition is not dirichlet"));
    Assert(dirichlet_value_map.contains(boundary_id),
           dealii::ExcMessage("Value does not exist in the map."));
    return dirichlet_value_map.at(boundary_id);
  }

  /**
   * \brief Add the value for a homogenous dirichlet boundary condition.
   */
  void
  add_boundary_condition(dealii::types::boundary_id boundary_id, double boundary_value)
  {
    Assert(boundary_condition_map.at(boundary_id) == type::DIRICHLET,
           dealii::ExcMessage("Boundary condition is not dirichlet"));
    dirichlet_value_map.emplace(boundary_id, boundary_value);
  }

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
        case type::TIME_DEPENDENT_NON_UNIFORM_DIRICHLET:
          return "TIME_DEPENDENT_NON_UNIFORM_DIRICHLET";
        case type::TIME_DEPENDENT_NON_UNIFORM_NEUMANN:
          return "TIME_DEPENDENT_NON_UNIFORM_NEUMANN";
        default:
          return "UNKNOWN";
      }
  }

private:
  // Map of boundary conditions and domain boundary for which they correspond to. For a
  // simple geometry like a square the boundary ids are marked, in order, by x=0, x=max,
  // y=0, y=max. More complex geometries can have somewhat arbitrary ordering, but will
  // render some of our assertions moot.
  std::map<dealii::types::boundary_id, type> boundary_condition_map;

  // A map of boundary values for dirichlet boundary conditions
  std::map<dealii::types::boundary_id, double> dirichlet_value_map;
};

/**
 * \brief Struct that holds boundary parameters.
 */
template <unsigned int dim>
struct BoundaryParameters
{
public:
  using BoundaryConditionMap =
    std::map<types::index, std::map<unsigned int, boundaryCondition>>;
  using BCList = std::map<types::index, std::map<unsigned int, std::string>>;
  using PinnedPointMap =
    std::map<types::index,
             std::pair<std::variant<double, std::vector<double>>, dealii::Point<dim>>>;

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

  /**
   * \brief Set the boundary condition string for a field index and component.
   */
  void
  set_boundary_condition_string(const std::string  &bc_string,
                                const types::index &index,
                                const unsigned int &component)
  {
    bc_list[index][component] = bc_string;
  }

  /**
   * \brief Whether there are time-dependent boundary conditions.
   */
  [[nodiscard]] bool
  has_time_dependent_bcs() const
  {
    return !time_dependent_bc_list.empty();
  }

  /**
   * \brief Whether the boundary condition is time-dependent.
   */
  [[nodiscard]] bool
  is_time_dependent(const types::index &index) const
  {
    return time_dependent_bc_list.contains(index);
  }

  /**
   * \brief Set a pinned point.
   */
  void
  set_pinned_point(const std::variant<double, std::vector<double>> &value,
                   const dealii::Point<dim>                        &point,
                   const types::index                              &index)
  {
    pinned_point_list[index] = std::make_pair(value, point);
  }

  /**
   * \brief Whether there is a pinned point for a field index.
   */
  [[nodiscard]] bool
  has_pinned_point(const types::index &index) const
  {
    return pinned_point_list.contains(index);
  }

  /**
   * \brief Get the pinned point for a field index.
   */
  [[nodiscard]] const std::pair<std::variant<double, std::vector<double>>,
                                dealii::Point<dim>> &
  get_pinned_point(const types::index &index) const
  {
    return pinned_point_list.at(index);
  }

  /**
   * \brief Get the boundary conditions list.
   */
  [[nodiscard]] const BoundaryConditionMap &
  get_boundary_condition_list() const
  {
    return boundary_condition_list;
  }

private:
  /**
   * \brief Set the boundary for a single component of a field index.
   */
  void
  set_boundary(const std::string  &bc_string,
               const types::index &index,
               const unsigned int &component);

  /**
   * \brief Perform a check on the boundary conditions to ensure that they are valid
   */
  void
  validate_boundary_conditions() const;

  // Map of unfiltered boundary conditions strings. The first key is the global index. The
  // second key is the number of dimensions.
  BCList bc_list;

  // Map of time-dependent boundary conditions strings. The first key is the global index.
  // The second key is the number of dimensions.
  std::set<types::index> time_dependent_bc_list;

  // Map of pinned points. The first key is the global index. The pair is the pinned value
  // and point.
  PinnedPointMap pinned_point_list = {};

  // Map of boundary conditions. The first key is the global index. The second key is the
  // number of dimensions.
  BoundaryConditionMap boundary_condition_list;
};

template <unsigned int dim>
inline void
BoundaryParameters<dim>::postprocess_and_validate(
  const std::map<unsigned int, variableAttributes> &var_attributes)
{
  for (const auto &[index, variable] : var_attributes)
    {
      // Ensure that boundary conditions are specified for all variables and their
      // components
      if (variable.get_field_type() == FieldType::VECTOR)
        {
          for (unsigned int i = 0; i < dim; i++)
            {
              if (variable.is_postprocess())
                {
                  set_boundary("NATURAL", variable.get_field_index(), i);
                }
              else
                {
                  AssertThrow(bc_list.contains(variable.get_field_index()),
                              dealii::ExcMessage("Invalid entry"));
                  AssertThrow(bc_list.at(variable.get_field_index()).contains(i),
                              dealii::ExcMessage("Invalid entry"));
                  AssertThrow(!bc_list.at(variable.get_field_index()).at(i).empty(),
                              dealii::ExcMessage(
                                "Boundary conditions must be specified "
                                "for all components in all vector field."));

                  set_boundary(bc_list.at(variable.get_field_index()).at(i),
                               variable.get_field_index(),
                               i);
                }
            }
        }
      else
        {
          if (variable.is_postprocess())
            {
              set_boundary("NATURAL", variable.get_field_index(), 0);
            }
          else
            {
              AssertThrow(bc_list.contains(variable.get_field_index()),
                          dealii::ExcMessage("Invalid entry"));
              AssertThrow(bc_list.at(variable.get_field_index()).contains(0),
                          dealii::ExcMessage("Invalid entry"));
              AssertThrow(!bc_list.at(variable.get_field_index()).at(0).empty(),
                          dealii::ExcMessage("Boundary conditions must be specified "
                                             "for all scalar fields."));

              set_boundary(bc_list.at(variable.get_field_index()).at(0),
                           variable.get_field_index(),
                           0);
            }
        }
    }

  // Validate boundary conditions
  validate_boundary_conditions();

  // Clear the bc_list now that it's no longer necessary
  bc_list.clear();

#ifdef ADDITIONAL_OPTIMIZATIONS
  // Check if any fields are duplicates in terms of boundary conditions
  // TODO (landinjm): Clean this up
  for (const auto &[index_1, variable_1] : var_attributes)
    {
      // Skip is the duplicate index has already been assigned
      if (variable_1.get_duplicate_field_index() != numbers::invalid_index)
        {
          continue;
        }
      if (variable_1.is_postprocess())
        {
          continue;
        }

      const auto field_type_1 = variable_1.get_field_type();

      for (const auto &[index_2, variable_2] : var_attributes)
        {
          if (variable_2.is_postprocess())
            {
              continue;
            }

          bool is_duplicate = false;

          const auto field_type_2 = variable_2.get_field_type();

          is_duplicate = field_type_1 == field_type_2 &&
                         check_duplicate_boundary_conditions(index_1, index_2);

          if (is_duplicate)
            {
              ConditionalOStreams::pout_verbose()
                << "Field " << variable_1.get_name()
                << " has the same boundary conditions as " << variable_2.get_name()
                << ". Using optimizations...\n";
              variable_2.set_duplicate_field_index(index_1);
            }
        }
    }
#endif
}

template <unsigned int dim>
inline bool
BoundaryParameters<dim>::check_duplicate_boundary_conditions(
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
BoundaryParameters<dim>::print_parameter_summary() const
{
  ConditionalOStreams::pout_summary()
    << "================================================\n"
    << "  Boundary Parameters\n"
    << "================================================\n";

  for (const auto &[index, component_map] : boundary_condition_list)
    {
      ConditionalOStreams::pout_summary() << "Index: " << index << "\n";
      for (const auto &[component, boundary_condition] : component_map)
        {
          ConditionalOStreams::pout_summary() << "  Component: " << component << "\n";
          for (const auto &[domain_id, boundary_type] :
               boundary_condition.get_boundary_condition_map())
            {
              ConditionalOStreams::pout_summary()
                << "    Boundary id: " << domain_id << "    "
                << "Condition: " << boundary_condition.to_string(boundary_type);
              if (boundary_type == boundaryCondition::type::DIRICHLET)
                {
                  ConditionalOStreams::pout_summary()
                    << " = " << boundary_condition.get_dirichlet_value(domain_id);
                }
              ConditionalOStreams::pout_summary() << "\n";
            }
        }
    }

  if (!pinned_point_list.empty())
    {
      ConditionalOStreams::pout_summary() << "Pinned field index: ";
    }
  for (const auto &[index, point_value_pair] : pinned_point_list)
    {
      ConditionalOStreams::pout_summary() << index << "\n"
                                          << "  Value: ";

      // Handle variant value printing
      std::visit(
        [&](const auto &value)
        {
          if constexpr (std::is_same_v<std::decay_t<decltype(value)>, double>)
            {
              ConditionalOStreams::pout_summary() << value;
            }
          else
            {
              for (unsigned int i = 0; i < value.size(); ++i)
                {
                  ConditionalOStreams::pout_summary() << value[i] << " ";
                }
            }
        },
        point_value_pair.first);

      ConditionalOStreams::pout_summary()
        << "\n  Point: " << point_value_pair.second << "\n";
    }
  ConditionalOStreams::pout_summary() << "\n" << std::flush;
}

template <unsigned int dim>
inline void
BoundaryParameters<dim>::set_boundary(const std::string  &bc_string,
                                      const types::index &index,
                                      const unsigned int &component)
{
  // Split string
  auto bc_string_list = dealii::Utilities::split_string_list(bc_string);

  // Check that there is either 1 or 2*dim entries in the vector. This can be changed
  // later to support other geometries.
  AssertThrow(bc_string_list.size() == 1 ||
                bc_string_list.size() == static_cast<std::size_t>(2 * dim),
              dealii::ExcMessage(
                "Either 1 or 2*dim boundary conditions must be specified."));

  // If there is only 1 boundary condition resize BC_string_list, copying the first
  // entry.
  if (bc_string_list.size() == 1)
    {
      bc_string_list.resize(static_cast<std::size_t>(2 * dim), bc_string_list[0]);
    }

  // Assign boundary condition
  boundaryCondition condition;
  for (unsigned int i = 0; i < (2 * dim); i++)
    {
      const std::string dirichlet = "DIRICHLET";
      const std::string neumann   = "NEUMANN";

      if (boost::iequals(bc_string_list[i], "NATURAL"))
        {
          condition.add_boundary_condition(i, boundaryCondition::type::NATURAL);
        }
      else if (boost::iequals(bc_string_list[i].substr(0, dirichlet.size()), dirichlet))
        {
          condition.add_boundary_condition(i, boundaryCondition::type::DIRICHLET);
          std::string dirichlet_value =
            bc_string_list[i].substr(dirichlet.size() + 1, bc_string_list[i].size());
          dirichlet_value = dealii::Utilities::trim(dirichlet_value);
          condition.add_boundary_condition(i,
                                           dealii::Utilities::string_to_double(
                                             dirichlet_value));
        }
      else if (boost::iequals(bc_string_list[i], "PERIODIC"))
        {
          condition.add_boundary_condition(i, boundaryCondition::type::PERIODIC);
        }
      else if (boost::iequals(bc_string_list[i].substr(0, neumann.size()), neumann))
        {
          AssertThrow(false, FeatureNotImplemented("Neumann boundary conditions"));
        }
      else if (boost::iequals(bc_string_list[i], "NON_UNIFORM_DIRICHLET"))
        {
          condition
            .add_boundary_condition(i, boundaryCondition::type::NON_UNIFORM_DIRICHLET);
        }
      else if (boost::iequals(bc_string_list[i], "NON_UNIFORM_NEUMANN"))
        {
          AssertThrow(false,
                      FeatureNotImplemented("Nonuniform neumann boundary conditions"));
        }
      else if (boost::iequals(bc_string_list[i], "TIME_DEPENDENT_NON_UNIFORM_DIRICHLET"))
        {
          condition.add_boundary_condition(
            i,
            boundaryCondition::type::TIME_DEPENDENT_NON_UNIFORM_DIRICHLET);
          time_dependent_bc_list.insert(index);
        }
      else if (boost::iequals(bc_string_list[i], "TIME_DEPENDENT_NON_UNIFORM_NEUMANN"))
        {
          AssertThrow(false,
                      FeatureNotImplemented(
                        "Time-dependent neumann boundary conditions"));
          time_dependent_bc_list.insert(index);
        }
      else
        {
          AssertThrow(false,
                      dealii::ExcMessage("Invalid boundary condition " +
                                         bc_string_list[i]));
        }
      // If periodic boundary conditions are used, ensure that they are applied on
      // both sides of the domain.
      if (i % 2 == 0)
        {
          AssertThrow(boost::iequals(bc_string_list[i], "PERIODIC") ==
                        boost::iequals(bc_string_list[i + 1], "PERIODIC"),
                      dealii::ExcMessage("Periodic boundary condition must be "
                                         "specified on both sides of domain"));
        }
    }

  boundary_condition_list[index].emplace(component, condition);
}

template <unsigned int dim>
inline void
BoundaryParameters<dim>::validate_boundary_conditions() const
{
  // Throw a warning if the pinned point is not on a vertex
  // TODO (landinjm): How do we want to handle this?
  for (const auto &[index, point_value_pair] : pinned_point_list)
    {
      Assert(point_value_pair.second == dealii::Point<dim>(),
             dealii::ExcMessage("Pinned point must be on the origin"));

      // Validate that vector values have the correct size
      std::visit(
        [&](const auto &value)
        {
          if constexpr (std::is_same_v<std::decay_t<decltype(value)>,
                                       std::vector<double>>)
            {
              AssertThrow(value.size() == dim,
                          dealii::ExcMessage("Vector value size must match dimension"));
            }
        },
        point_value_pair.first);
    }

  // Throw a warning if only some fields have periodic boundary conditions
  std::vector<bool> periodic_ids(static_cast<dealii::types::boundary_id>(2 * dim), false);
  for (const auto &[index, component_map] : boundary_condition_list)
    {
      for (const auto &[component, boundary_condition] : component_map)
        {
          for (const auto &[domain_id, boundary_type] :
               boundary_condition.get_boundary_condition_map())
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
               boundary_condition.get_boundary_condition_map())
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