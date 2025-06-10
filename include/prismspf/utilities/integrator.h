// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Compute the integral of a given field.
 */
template <unsigned int dim, unsigned int degree, typename number>
class Integrator
{
public:
  using VectorType = dealii::LinearAlgebra::distributed::Vector<number>;

  /**
   * \brief Constructor.
   */
  Integrator() = default;

  /**
   * \brief Compute the integral for a scalar field
   */
  void
  compute_integral(number                        &integral_value,
                   const dealii::DoFHandler<dim> &dof_handler,
                   const VectorType              &vector) const;

  /**
   * \brief Compute the integral for a vector field
   */
  void
  compute_integral(std::vector<number>           &integral_value,
                   const dealii::DoFHandler<dim> &dof_handler,
                   const VectorType              &vector) const;
};

PRISMS_PF_END_NAMESPACE
