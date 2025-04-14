// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
class userInputParameters;

/**
 * \brief Forward declaration of user-facing implementation
 */
template <int dim>
class customNonuniformDirichlet;

/**
 * \brief Function for user-implemented nonuniform dirichlet boundary condition.
 */
template <int dim, fieldType field_type = fieldType::SCALAR>
class nonuniformDirichlet : public dealii::Function<dim, double>
{
public:
  /**
   * \brief Constructor.
   */
  nonuniformDirichlet(const unsigned int             &_index,
                      const unsigned int             &_boundary_id,
                      const userInputParameters<dim> &_user_inputs);

  // NOLINTBEGIN(readability-identifier-length, readability-avoid-const-params-in-decls)

  /**
   * \brief Scalar value.
   */
  double
  value(const dealii::Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * \brief Vector value.
   */
  void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<double> &value) const override;

  // NOLINTEND(readability-identifier-length, readability-avoid-const-params-in-decls)

private:
  unsigned int index;

  unsigned int boundary_id;

  const userInputParameters<dim> *user_inputs;

  customNonuniformDirichlet<dim> custom_nonuniform_dirichlet;
};

/**
 * \brief User-facing implementation of nonuniform boundary conditions
 */
template <int dim>
class customNonuniformDirichlet
{
public:
  /**
   * \brief Constructor.
   */
  customNonuniformDirichlet() = default;

  /**
   * \brief Function that passes the value/vector and point that are set in the nonuniform
   * dirichlet.
   */
  void
  set_nonuniform_dirichlet(const unsigned int             &index,
                           const unsigned int             &boundary_id,
                           const unsigned int             &component,
                           const dealii::Point<dim>       &point,
                           double                         &scalar_value,
                           double                         &vector_component_value,
                           const userInputParameters<dim> &user_inputs) const;
};

PRISMS_PF_END_NAMESPACE
