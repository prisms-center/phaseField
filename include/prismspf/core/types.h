// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/config.h>

#include <limits>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * A collection of types that we use throughout the code.
 */
namespace Types
{
  /**
   * @brief Type for field indices.
   */
  using Index = unsigned int;

} // namespace Types

/**
 * A collection of numbers we use throughout the code.
 */
namespace Numbers
{
  /**
   * @brief Invalid field index.
   */
  static constexpr Types::Index invalid_index = -1;

  /**
   * @brief Max number of subsections.
   */
  static constexpr unsigned int default_subsections = 5;

} // namespace Numbers

/**
 * A collection of defaults we use throughout the code.
 */
namespace Defaults
{
  /**
   * @brief Machine epsilon
   *
   * @todo Add a description of when and why we use this
   */
  template <typename RealType>
  static constexpr RealType machine_epsilon = std::numeric_limits<RealType>::epsilon();

  /**
   * @brief Tolerance factor
   *
   * @todo Add a description of when and why we use this
   */
  template <typename RealType>
  static constexpr RealType tolerance = RealType(10) * machine_epsilon<RealType>;

} // namespace Defaults

/**
 * A collection of global type definitions that we use elsewhere in the code.
 */

using EvalFlags = dealii::EvaluationFlags::EvaluationFlags;

PRISMS_PF_END_NAMESPACE
