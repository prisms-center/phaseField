// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/types.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Collection of system-wide shared objects.
 *
 * Provides access to objects that are used throughout the application and are expensive
 * to construct or identical for all users of a given template specialization. These
 * objects are stored once per `(dim, degree)` specialization and reused wherever needed
 * to avoid repeated construction and ensure consistent configuration.
 *
 * Typical examples include quadrature rules, mappings, finite elements, and
 * other immutable data structures shared across the code base.
 *
 * We also use lazy initializers for the expensive objects that aren't always
 * constructors.
 */
template <unsigned int dim, unsigned int degree>
class SystemWide
{
public:
  /**
   * @brief Scalar and vector FE systems
   */
  inline static const std::array<const dealii::FESystem<dim>, 2> fe_systems = {
    dealii::FESystem<dim>(dealii::FE_Q<dim>(dealii::QGaussLobatto<1>(degree + 1)), 1),
    dealii::FESystem<dim>(dealii::FE_Q<dim>(dealii::QGaussLobatto<1>(degree + 1)), dim)};

  /**
   * @brief Mappings to and from reference cell
   */
  inline static const dealii::MappingQ1<dim> mapping;

  /**
   * @brief Quadrature rule
   */
  inline static const dealii::QGaussLobatto<dim> quadrature =
    dealii::QGaussLobatto<dim>(degree + 1);

  /**
   * @brief Quadrature rule for MatrixFree
   *
   * @todo Figure out why dim = 1
   */
  inline static const dealii::QGaussLobatto<1> quadrature_matrix_free =
    dealii::QGaussLobatto<1>(degree + 1);
};

PRISMS_PF_END_NAMESPACE
