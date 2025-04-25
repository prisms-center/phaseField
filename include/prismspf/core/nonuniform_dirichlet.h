// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class userInputParameters;

/**
 * \brief Forward declaration of user-facing implementation
 */
template <unsigned int dim, typename number>
class customNonuniformDirichlet;

/**
 * \brief Function for user-implemented nonuniform dirichlet boundary condition.
 */
template <unsigned int dim, typename number>
class nonuniformDirichlet : public dealii::Function<dim, number>
{
public:
  /**
   * \brief Constructor.
   */
  nonuniformDirichlet(unsigned int                    _index,
                      unsigned int                    _boundary_id,
                      const userInputParameters<dim> &_user_inputs,
                      unsigned int                    spacedim);

  // NOLINTBEGIN(readability-identifier-length, readability-avoid-const-params-in-decls)

  /**
   * \brief Scalar value.
   */
  number
  value(const dealii::Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * \brief Vector value.
   */
  void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<number> &value) const override;

  // NOLINTEND(readability-identifier-length, readability-avoid-const-params-in-decls)

private:
  unsigned int index;

  unsigned int boundary_id;

  const userInputParameters<dim> *user_inputs;

  customNonuniformDirichlet<dim, number> custom_nonuniform_dirichlet;
};

/**
 * \brief User-facing implementation of nonuniform boundary conditions
 */
template <unsigned int dim, typename number>
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
                           number                         &scalar_value,
                           number                         &vector_component_value,
                           const userInputParameters<dim> &user_inputs) const;
};

PRISMS_PF_END_NAMESPACE
