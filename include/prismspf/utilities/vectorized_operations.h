// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/vectorization.h>

#include <prismspf/config.h>

/**
 * This file provides some operations on dealii::VectorizedArray that are not provided by
 * dealii
 */
namespace dealii
{

  template <typename number, typename other_number>
  dealii::VectorizedArray<number>
  fmod(const dealii::VectorizedArray<number> &value, const other_number &modulus)
  {
    using std::fmod;
    dealii::VectorizedArray<number> out;
    for (unsigned int index = 0; index < dealii::VectorizedArray<number>::size(); ++index)
      {
        out[index] = fmod(value[index], modulus);
      }
    return out;
  }

  template <typename number, typename other_number>
  dealii::VectorizedArray<number>
  fmod(const dealii::VectorizedArray<number>       &value,
       const dealii::VectorizedArray<other_number> &modulus)
  {
    using std::fmod;
    dealii::VectorizedArray<number> out;
    for (unsigned int index = 0; index < dealii::VectorizedArray<number>::size(); ++index)
      {
        out[index] = fmod(value[index], modulus[index]);
      }
    return out;
  }
} // namespace dealii