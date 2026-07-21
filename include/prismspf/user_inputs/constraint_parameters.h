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
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/user_inputs/parameter_base.h>

#include <prismspf/config.h>

#include <map>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Condition of boundary condition.
 */
enum Condition : std::uint8_t
{
  Natural,
  Dirichlet,
  Neumann,
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
  if (boundary_string == "Periodic")
    {
      return Condition::Periodic;
    }
  AssertThrow(false, dealii::ExcMessage("Invalid boundary condition " + boundary_string));
  return Condition::Natural;
}

/**
 * @brief Map of boundary conditions and domain boundary for which they correspond to.
 * For a simple geometry like a square the boundary ids are marked, in order, by x=0,
 * x=max, y=0, y=max.
 */
using ComponentConditions = std::map<unsigned int, Condition>;

struct BoundaryConditionSet
{
  std::map<unsigned int, ComponentConditions> component_constraints;

  bool time_dependent = false;
};

/**
 * @brief Struct that holds boundary parameters.
 */
struct BoundaryParameters
{
  /**
   * @brief Declare the parameters to be read from file.
   */
  static void
  declare(dealii::ParameterHandler &parameter_handler,
          unsigned int              n_subsections = Numbers::default_subsections);

  /**
   * @brief Assign the parameters from file.
   */
  template <unsigned int dim = 3>
  void
  assign(dealii::ParameterHandler &parameter_handler,
         unsigned int              n_subsections = Numbers::default_subsections);

  /**
   * @brief Validate.
   */
  void
  validate(const std::vector<FieldAttributes> &field_attributes,
           const std::vector<SolveBlock>      &solve_blocks) const;

  // Map of boundary conditions. The first key is the field name.
  std::unordered_map<std::string, BoundaryConditionSet> boundary_condition_list;
};

PRISMS_PF_END_NAMESPACE
