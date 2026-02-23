// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/point.h>
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
#include <cstddef>
#include <map>
#include <string>
#include <variant>
#include <vector>

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
[[nodiscard]] std::string
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
[[nodiscard]] Condition
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
public:
  /**
   * @brief Test for equality of two boundary conditions.
   */
  bool
  operator==(const ComponentConditions &other) const
  {
    return conditions == other.conditions &&
           uniform_dirichlet_values == other.uniform_dirichlet_values;
  }

  /**
   * @brief Get the map of boundary conditions.
   */
  [[nodiscard]] const std::vector<Condition> &
  get_conditions() const
  {
    return conditions;
  }

  /**
   * @brief Get the map of boundary conditions.
   */
  [[nodiscard]] Condition
  get_condition(unsigned int boundary_id) const
  {
    return conditions[boundary_id];
  }

  /**
   * @brief Add a boundary conditions.
   */
  void
  add_boundary_condition(unsigned int boundary_id, Condition boundary_type)
  {
    conditions[boundary_id] = boundary_type;
  }

  /**
   * @brief Get the value for a homogeneous dirichlet boundary condition.
   */
  [[nodiscard]] double
  get_uniform_dirichlet_value(unsigned int boundary_id) const
  {
    Assert(uniform_dirichlet_values.size() > boundary_id,
           dealii::ExcMessage("Boundary condition is not dirichlet"));
    return uniform_dirichlet_values.at(boundary_id);
  }

  /**
   * @brief Add the value for a homogeneous dirichlet boundary condition.
   */
  void
  add_uniform_dirichlet_value(unsigned int boundary_id, double boundary_value)
  {
    uniform_dirichlet_values[boundary_id] = boundary_value;
  }

private:
  // Map of boundary conditions and domain boundary for which they correspond to. For a
  // simple geometry like a square the boundary ids are marked, in order, by x=0, x=max,
  // y=0, y=max. More complex geometries can have somewhat arbitrary ordering, but will
  // render some of our assertions moot.
  std::map<unsigned int, Condition> conditions;

  // A map of boundary values for dirichlet boundary conditions
  std::map<unsigned int, double> uniform_dirichlet_values;
};

/**
 * @brief Struct that holds boundary parameters.
 */
template <unsigned int dim>
struct BoundaryParameters
{
public:
  using FieldConstraints = std::array<ComponentConditions, dim>;
  using PinnedPointMap   = std::map<unsigned int, dealii::Point<dim>>;

  /**
   * @brief Postprocess and validate parameters.
   */
  void
  postprocess_and_validate(std::vector<FieldAttributes> &field_attributes);

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
  {
    return !time_dependent_bc_list.empty();
  }

  /**
   * @brief Set a pinned point.
   */
  void
  set_pinned_point(const std::variant<double, std::vector<double>> &value,
                   const dealii::Point<dim>                        &point,
                   const Types::Index                              &index)
  {
    pinned_point_list[index] = std::make_pair(value, point);
  }

  /**
   * @brief Whether there is a pinned point for a field index.
   */
  [[nodiscard]] bool
  has_pinned_point(const Types::Index &index) const
  {
    return pinned_point_list.contains(index);
  }

  /**
   * @brief Get the pinned point for a field index.
   */
  [[nodiscard]] const std::pair<std::variant<double, std::vector<double>>,
                                dealii::Point<dim>> &
  get_pinned_point(const Types::Index &index) const
  {
    return pinned_point_list.at(index);
  }

  /**
   * @brief Get the boundary conditions list.
   */
  [[nodiscard]] const std::map<unsigned int, FieldConstraints> &
  get_boundary_conditions() const
  {
    return boundary_condition_list;
  }

  /**
   * @brief Check if any boundaries are periodic
   */
  [[nodiscard]] bool
  has_periodic_boundaries() const
  {
    return std::any_of(boundary_condition_list.begin(),
                       boundary_condition_list.end(),
                       [](const FieldConstraints &val)
                       {

                       });
  }

private:
  /**
   * @brief Set the boundary for a single component of a field index.
   */
  void
  set_boundary(const Types::Index &index,
               const unsigned int &component,
               const std::string  &bc_string);

  /**
   * @brief Perform a check on the boundary conditions to ensure that they are valid
   */
  void
  validate_boundary_conditions() const;

  // Map of pinned points. The first key is the global index. The pair is the pinned
  // value and point.
  PinnedPointMap pinned_point_list = {};

  // Map of boundary conditions. The first key is the field index. The std::vector index
  // is the component of the field ({0} for scalars, {1,2,3} for vectors)
  std::map<unsigned int, FieldConstraints> boundary_condition_list;
};

template <unsigned int dim>
inline void
BoundaryParameters<dim>::postprocess_and_validate(
  std::vector<FieldAttributes> &field_attributes)
{}

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
              if (boundary_type == Condition::Dirichlet)
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
                                      const Types::Index &index,
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
  BoundaryCondition condition;
  for (unsigned int i = 0; i < (2 * dim); i++)
    {
      const std::string dirichlet = "Dirichlet";
      const std::string neumann   = "Neumann";

      if (boost::iequals(bc_string_list[i], "Natural"))
        {
          condition.add_boundary_condition(i, Condition::Natural);
        }
      else if (boost::iequals(bc_string_list[i].substr(0, dirichlet.size()), dirichlet))
        {
          condition.add_boundary_condition(i, Condition::Dirichlet);
          std::string dirichlet_value =
            bc_string_list[i].substr(dirichlet.size() + 1, bc_string_list[i].size());
          dirichlet_value = dealii::Utilities::trim(dirichlet_value);
          condition.add_boundary_condition(i,
                                           dealii::Utilities::string_to_double(
                                             dirichlet_value));
        }
      else if (boost::iequals(bc_string_list[i], "Periodic"))
        {
          condition.add_boundary_condition(i, Condition::Periodic);
        }
      else if (boost::iequals(bc_string_list[i].substr(0, neumann.size()), neumann))
        {
          AssertThrow(false, FeatureNotImplemented("Neumann boundary conditions"));
        }
      else if (boost::iequals(bc_string_list[i], "NonuniformDirichlet"))
        {
          condition.add_boundary_condition(i, Condition::NonuniformDirichlet);
        }
      else if (boost::iequals(bc_string_list[i], "NonuniformNeumann"))
        {
          AssertThrow(false,
                      FeatureNotImplemented("Nonuniform neumann boundary conditions"));
        }
      else if (boost::iequals(bc_string_list[i], "TimeDependentNonuniformDirichlet"))
        {
          condition.add_boundary_condition(i,
                                           Condition::TimeDependentNonuniformDirichlet);
          time_dependent_bc_list.insert(index);
        }
      else if (boost::iequals(bc_string_list[i], "TimeDependentNonuniformNeumann"))
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
          AssertThrow(boost::iequals(bc_string_list[i], "Periodic") ==
                        boost::iequals(bc_string_list[i + 1], "Periodic"),
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
      // Disable this TODO (landinjm): Fix
      /*Assert(point_value_pair.second == dealii::Point<dim>(),
             dealii::ExcMessage("Pinned point must be on the origin"));*/

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
              if (boundary_type == Condition::Periodic)
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
              if (boundary_type != Condition::Periodic && periodic_ids[domain_id])
                {
                  // TODO (landinjm): This needs to be fixed to ignore postprocess
                  // fields. Disable it for now.
                  /*AssertThrow(
                    false,
                    dealii::ExcMessage(
                      "All fields for a given domain id (side) must have periodic "
                      "boundary conditions if any field has periodic boundary "
                      "conditions"));*/
                }
            }
        }
    }
}

PRISMS_PF_END_NAMESPACE
