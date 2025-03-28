// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

namespace types
{
  /**
   * \brief Type for field indices.
   */
  using index = unsigned int;

} // namespace types

namespace numbers
{
  /**
   * \brief Invalid field index.
   */
  static const types::index invalid_index = static_cast<types::index>(-1);

} // namespace numbers

namespace defaults
{
  /**
   * \brief Default field index.
   */
  static const types::index index = 0;

  /**
   * \brief Default tolerance.
   */
  static const double tolerance = 1.0e-6;

  /**
   * \brief Default iterations.
   */
  static const unsigned int iterations = 100;

  /**
   * \brief Default eigenvalue smoothing range for multigrid.
   */
  static const double smoothing_range = 15.0;

  /**
   * \brief Default smoother degree for multigrid.
   */
  static const unsigned int smoother_degree = 5;

  /**
   * \brief Default CG iterations to find the maximum eigenvalue for multigrid.
   */
  static const unsigned int eig_cg_n_iterations = 10;

} // namespace defaults

PRISMS_PF_END_NAMESPACE