// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/matrix_free/evaluation_flags.h>

#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

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
   * @brief Max number of subsections.
   */
  static constexpr unsigned int default_subsections = 5;

} // namespace Numbers

namespace Defaults
{} // namespace Defaults

using EvalFlags = dealii::EvaluationFlags::EvaluationFlags;

PRISMS_PF_END_NAMESPACE
