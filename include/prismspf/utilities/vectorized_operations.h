// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/vectorization.h>

#include <prismspf/config.h>

#include <cmath>

/**
 * This file provides some operations on dealii::VectorizedArray that are not provided by
 * dealii. Use the std namespace so we can call them with std::function
 */
namespace std
{
  // NOLINTBEGIN(cert-dcl58-cpp, readability-identifier-naming,
  // readability-identifier-length)
  // clang-format off

  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  erf(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size(); ++i)
      out[i] = std::erf(x[i]);
    return out;
  }

  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  erfc(const ::dealii::VectorizedArray<Number, width> &x)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size(); ++i)
      out[i] = std::erfc(x[i]);
    return out;
  }

  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  fmod(const ::dealii::VectorizedArray<Number, width> &numer, const Number denom)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size(); ++i)
      out[i] = std::fmod(numer[i], denom);
    return out;
  }

  template <typename Number, std::size_t width>
  inline ::dealii::VectorizedArray<Number, width>
  fmod(const ::dealii::VectorizedArray<Number, width> &numer,
       const ::dealii::VectorizedArray<Number, width> &denom)
  {
    ::dealii::VectorizedArray<Number, width> out;
    for (unsigned int i = 0; i < dealii::VectorizedArray<Number, width>::size(); ++i)
      out[i] = std::fmod(numer[i], denom[i]);
    return out;
  }

  // clang-format on
  // NOLINTEND(cert-dcl58-cpp, readability-identifier-naming,
  // readability-identifier-length)

} // namespace std
