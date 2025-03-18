// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef types_h
#define types_h

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

} // namespace defaults

PRISMS_PF_END_NAMESPACE

#endif