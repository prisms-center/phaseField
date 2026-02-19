// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/types.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This is the main class that handles the construction and solving of
 * user-specified PDEs.
 */
template <unsigned int dim, unsigned int degree>
class SystemWide
{
public:
  /**
   * @brief Scalar and Vector FE systems.
   */
  inline static const std::array<const dealii::FESystem<dim>, 2> fe_systems = {
    dealii::FESystem<dim>(dealii::FE_Q<dim>(dealii::QGaussLobatto<1>(degree + 1)), 1),
    dealii::FESystem<dim>(dealii::FE_Q<dim>(dealii::QGaussLobatto<1>(degree + 1)), dim)};

  /**
   * @brief Mappings to and from reference cell.
   */
  inline static const dealii::MappingQ1<dim> mapping;

  /**
   * @brief Quadrature rule
   */
  inline static const dealii::QGaussLobatto<dim> quadrature =
    dealii::QGaussLobatto<dim>(degree + 1);
};

PRISMS_PF_END_NAMESPACE