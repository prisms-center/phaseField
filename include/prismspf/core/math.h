#pragma once

#include <prismspf/config.h>

#ifdef PRISMS_PF_GPU
#  include <Kokkos_Core.hpp>
#  include <Kokkos_MathematicalFunctions.hpp>
#else
#  include <cmath>
#endif

PRISMS_PF_BEGIN_NAMESPACE

namespace math
{
#ifdef PRISMS_PF_GPU
#  define PRISMS_PF_MATHTYPE Kokkos
#else
#  define PRISMS_PF_MATHTYPE std
#endif

  // Basic operations
  using PRISMS_PF_MATHTYPE::fabs;
  using PRISMS_PF_MATHTYPE::fdim;
  using PRISMS_PF_MATHTYPE::fmax;
  using PRISMS_PF_MATHTYPE::fmin;
  using PRISMS_PF_MATHTYPE::fmod;
  using PRISMS_PF_MATHTYPE::remainder;

  // Exponentials and logarithms
  using PRISMS_PF_MATHTYPE::exp;
  using PRISMS_PF_MATHTYPE::log;
  using PRISMS_PF_MATHTYPE::log10;
  using PRISMS_PF_MATHTYPE::log1p;
  using PRISMS_PF_MATHTYPE::log2;

  // Power functions
  using PRISMS_PF_MATHTYPE::cbrt;
  using PRISMS_PF_MATHTYPE::pow;
  using PRISMS_PF_MATHTYPE::sqrt;

  // Trigonometric functions
  using PRISMS_PF_MATHTYPE::acos;
  using PRISMS_PF_MATHTYPE::asin;
  using PRISMS_PF_MATHTYPE::atan;
  using PRISMS_PF_MATHTYPE::atan2;
  using PRISMS_PF_MATHTYPE::cos;
  using PRISMS_PF_MATHTYPE::sin;
  using PRISMS_PF_MATHTYPE::tan;

  // Hyperbolic functions
  using PRISMS_PF_MATHTYPE::acosh;
  using PRISMS_PF_MATHTYPE::asinh;
  using PRISMS_PF_MATHTYPE::atanh;
  using PRISMS_PF_MATHTYPE::cosh;
  using PRISMS_PF_MATHTYPE::sinh;
  using PRISMS_PF_MATHTYPE::tanh;

  // Error and gamma functions
  using PRISMS_PF_MATHTYPE::erf;
  using PRISMS_PF_MATHTYPE::erfc;
  using PRISMS_PF_MATHTYPE::lgamma;
  using PRISMS_PF_MATHTYPE::tgamma;

  // Nearest integer floating point operations
  using PRISMS_PF_MATHTYPE::ceil;
  using PRISMS_PF_MATHTYPE::floor;
  using PRISMS_PF_MATHTYPE::round;
  using PRISMS_PF_MATHTYPE::trunc;

#undef PRISMS_PF_MATHTYPE

} // namespace math

PRISMS_PF_END_NAMESPACE