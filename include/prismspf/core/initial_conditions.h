// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class userInputParameters;

/**
 * \brief Function for user-implemented initial conditions. These are only ever calculated
 * for explicit time dependent fields and implicit time dependent, as all others are
 * calculated at runtime.
 */
template <unsigned int dim, unsigned int degree>
class initialCondition : public dealii::Function<dim, double>
{
public:
  /**
   * \brief Constructor.
   */
  initialCondition(
    const unsigned int                                            &_index,
    const fieldType                                               &field_type,
    const std::shared_ptr<const PDEOperator<dim, degree, double>> &_pde_operator);

  // NOLINTBEGIN(readability-identifier-length)

  /**
   * \brief Scalar/Vector value.
   */
  void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<double> &value) const override;

  // NOLINTEND(readability-identifier-length)

private:
  unsigned int index;

  std::shared_ptr<const PDEOperator<dim, degree, double>> pde_operator;
};

///**
// * \brief User-facing implementation of initial conditions
// */
// template <unsigned int dim>
// class customInitialCondition
//{
// public:
//  /**
//   * \brief Constructor.
//   */
//  customInitialCondition() = default;
//
//  /**
//   * \brief Function that passes the value/vector and point that are set in the initial
//   * condition.
//   */
//  void
//  set_initial_condition(const unsigned int             &index,
//                        const unsigned int             &component,
//                        const dealii::Point<dim>       &point,
//                        double                         &scalar_value,
//                        double                         &vector_component_value,
//                        const userInputParameters<dim> &user_inputs) const;
//};

PRISMS_PF_END_NAMESPACE
