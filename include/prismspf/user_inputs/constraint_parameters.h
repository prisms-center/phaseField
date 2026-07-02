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

#include <prismspf/user_inputs/parameter_base.h>

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
  has_time_dependent_bcs() const;
};

template <unsigned int dim>
struct FieldConstraints
{
  std::array<ComponentConditions, dim> component_constraints;

  [[nodiscard]] bool
  has_time_dependent_bcs() const;
};

/**
 * @brief Struct that holds boundary parameters.
 */
template <unsigned int dim>
struct BoundaryParameters : public ParameterBase
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  predeclare(dealii::ParameterHandler &parameter_handler) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  preassign(dealii::ParameterHandler &parameter_handler) override;

  /**
   * @brief Declare the parameters to be read from file.
   */
  void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int max_criteria = Numbers::max_subsections) const override;

  /**
   * @brief Assign the parameters from file.
   */
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              max_criteria = Numbers::max_subsections) override;

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const override;

  // Map of boundary conditions. The first key is the field index.
  std::unordered_map<std::string, FieldConstraints<dim>> boundary_condition_list;
};

PRISMS_PF_END_NAMESPACE
