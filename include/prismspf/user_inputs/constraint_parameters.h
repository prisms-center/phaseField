// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

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
#include <set>
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
  // Map of boundary conditions and domain boundary for which they correspond to. For a
  // simple geometry like a square the boundary ids are marked, in order, by x=0, x=max,
  // y=0, y=max. More complex geometries can have somewhat arbitrary ordering, but will
  // render some of our assertions moot.
  std::map<unsigned int, Condition> conditions;
};

template <unsigned int dim>
struct FieldConstraints
{
  std::array<ComponentConditions, dim>                            component_constraints;
  std::set<std::pair<dealii::Point<dim>, dealii::Tensor<1, dim>>> pinned_points;
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
  validate(const std::vector<FieldAttributes> &field_attributes);

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
  }

  // Map of boundary conditions. The first key is the field index.
  std::map<std::string, FieldConstraints<dim>> boundary_condition_list;
};

template <unsigned int dim>
inline void
BoundaryParameters<dim>::validate(
  [[maybe_unused]] const std::vector<FieldAttributes> &field_attributes)
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

PRISMS_PF_END_NAMESPACE
