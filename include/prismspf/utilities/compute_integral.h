// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/vector.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Compute the integral of a given field.
 */
template <int dim, int degree, typename number>
class computeIntegral
{
public:
  using VectorType = dealii::LinearAlgebra::distributed::Vector<number>;

  /**
   * \brief Constructor.
   */
  explicit computeIntegral() = default;

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
