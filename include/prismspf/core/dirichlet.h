// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/function.h>

#include <prismspf/core/simulation_timer.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

template <unsigned int dim, unsigned int degree, typename number>
class PDEOperatorBase;

/**
 * @brief Function for user-implemented Dirichlet boundary condition.
 */
template <unsigned int dim, unsigned int degree, typename number>
class DirichletConditions : public dealii::Function<dim, number>
{
public:
  /**
   * @brief Constructor.
   */
  DirichletConditions(unsigned int                                _index,
                      unsigned int                                _boundary_id,
                      const PDEOperatorBase<dim, degree, number> &_pde_operator,
                      const SimulationTimer                      &_sim_timer,
                      unsigned int                                spacedim);

  // NOLINTBEGIN(readability-identifier-length, readability-avoid-const-params-in-decls)

  /**
   * @brief Scalar value.
   */
  number
  value(const dealii::Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * @brief Vector value.
   */
  void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<number> &value) const override;

  // NOLINTEND(readability-identifier-length, readability-avoid-const-params-in-decls)

private:
  unsigned int index;

  unsigned int boundary_id;

  const PDEOperatorBase<dim, degree, number> *pde_operator = nullptr;

  const SimulationTimer *sim_timer = nullptr;
};

PRISMS_PF_END_NAMESPACE
