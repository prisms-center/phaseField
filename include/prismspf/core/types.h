// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

namespace Types
{
  /**
   * @brief Type for field indices.
   */
  using Index = unsigned int;

} // namespace Types

namespace Numbers
{
  /**
   * @brief Invalid field index.
   */
  static const Types::Index invalid_index = static_cast<Types::Index>(-1);

  /**
   * @brief Max element degree.
   */
  static const unsigned int max_element_degree = 6;

  /**
   * @brief Max number of subsections.
   */
  static const unsigned int max_subsections = 16;

  /**
   * @brief Invalid PDE type.
   */
  static const PDEType invalid_pde_type = static_cast<PDEType>(-1);

  /**
   * @brief Invalid field solve type.
   */
  static const FieldSolveType invalid_field_solve_type = static_cast<FieldSolveType>(-1);

} // namespace Numbers

namespace Defaults
{
  /**
   * @brief Default field index.
   */
  static const Types::Index index = 0;

  /**
   * @brief Default tolerance.
   */
  static const double tolerance = 1.0e-6;

  /**
   * @brief Default mesh tolerance.
   */
  static const double mesh_tolerance = 1.0e-15;

  /**
   * @brief Default iterations.
   */
  static const unsigned int iterations = 100;

  /**
   * @brief Default eigenvalue smoothing range for multigrid.
   */
  static const double smoothing_range = 15.0;

  /**
   * @brief Default smoother degree for multigrid.
   */
  static const unsigned int smoother_degree = 5;

  /**
   * @brief Default CG iterations to find the maximum eigenvalue for multigrid.
   */
  static const unsigned int eig_cg_n_iterations = 10;

} // namespace Defaults

PRISMS_PF_END_NAMESPACE
