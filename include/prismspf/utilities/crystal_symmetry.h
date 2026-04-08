#pragma once

#include <deal.II/base/vectorization.h>

#include <prismspf/utilities/vectorized_operations.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

namespace Symmetries
{
  /**
   * @brief Compute cos(N * arctan(x)).
   *
   * This function provides an efficient way to calculate symmetries like $cos(N
   * arctan(\theta))$, where `N` is a known number at compile time.
   *
   * Importantly, it doesn't always make sense to unroll the function with a Chebyshev
   * polynomial due to the crossover point in multiplication operations and using the math
   * library. We've found this to occur at N = 3, although your mileage will vary based on
   * your architecture.
   * TODO: Benchmark this better
   */
  template <unsigned int N, typename T>
  [[nodiscard]] inline DEAL_II_ALWAYS_INLINE T
  cos_arctan(const T &x) noexcept
  {
    using std::atan;
    using std::cos;
    using std::sqrt;

    if constexpr (N <= 3)
      {
        const T t1 = T(1) / sqrt(T(1) + x * x);
        const T t2 = t1 * t1;

        if constexpr (N == 0)
          {
            return T(1);
          }
        else if constexpr (N == 1)
          {
            return t1;
          }
        else if constexpr (N == 2)
          {
            return T(2) * t2 - T(1);
          }
        else if constexpr (N == 3)
          {
            return t1 * (T(4) * t2 - T(3));
          }
        else if constexpr (N == 4)
          {
            return T(8) * t2 * t2 - T(8) * t2 + T(1);
          }
      }
    else
      {
        return std::cos(T(N) * std::atan(x));
      }
  }

  /**
   * @brief Compute cos(N * theta) where theta = arctan(ny/nx).
   *
   * This function provides an efficient way to calculate symmetries like $cos(N
   * arctan(n_y/n_x))$, where `N` is a known number at compile time.
   *
   * Importantly, it doesn't always make sense to unroll the function with a
   * Chebyshev polynomial due to the crossover point in multiplication operations and
   * using the math library. We've found this to occur at N = 5, although your mileage
   * will vary based on your architecture.
   * TODO: Benchmark this better
   */
  template <unsigned int N, typename T>
  [[nodiscard]] inline DEAL_II_ALWAYS_INLINE T
  cos_theta(const T &nx, const T &ny) noexcept
  {
    using std::atan2;
    using std::cos;

    if constexpr (N <= 5)
      {
        const T t1  = nx;      // cos θ
        const T t2  = t1 * t1; // (cos θ)^2
        const T t3  = t2 * t1; // (cos θ)^3
        const T t4  = t2 * t2; // (cos θ)^4
        const T t5  = t3 * t2; // (cos θ)^5
        const T t6  = t3 * t3; // (cos θ)^6
        const T t7  = t4 * t3; // (cos θ)^7
        const T t8  = t4 * t4; // (cos θ)^8
        const T t9  = t5 * t4; // (cos θ)^9
        const T t10 = t5 * t5; // (cos θ)^10
        const T t11 = t6 * t5; // (cos θ)^11
        const T t12 = t6 * t6; // (cos θ)^12

        if constexpr (N == 0)
          {
            return T(1);
          }
        else if constexpr (N == 1)
          {
            return t1;
          }
        else if constexpr (N == 2)
          {
            return T(2) * t2 - T(1);
          }
        else if constexpr (N == 3)
          {
            return T(4) * t3 - T(3) * t1;
          }
        else if constexpr (N == 4)
          {
            return T(8) * t4 - T(8) * t2 + T(1);
          }
        else if constexpr (N == 5)
          {
            return T(16) * t5 - T(20) * t3 + T(5) * t1;
          }
        else if constexpr (N == 6)
          {
            return T(32) * t6 - T(48) * t4 + T(18) * t2 - T(1);
          }
        else if constexpr (N == 7)
          {
            return T(64) * t7 - T(112) * t5 + T(56) * t3 - T(7) * t1;
          }
        else if constexpr (N == 8)
          {
            return T(128) * t8 - T(256) * t6 + T(160) * t4 - T(32) * t2 + T(1);
          }
        else if constexpr (N == 9)
          {
            return T(256) * t9 - T(576) * t7 + T(432) * t5 - T(120) * t3 + T(9) * t1;
          }
        else if constexpr (N == 10)
          {
            return T(512) * t10 - T(1280) * t8 + T(1120) * t6 - T(400) * t4 + T(50) * t2 -
                   T(1);
          }
        else if constexpr (N == 11)
          {
            return T(1024) * t11 - T(2816) * t9 + T(2816) * t7 - T(1232) * t5 +
                   T(220) * t3 - T(11) * t1;
          }
        else if constexpr (N == 12)
          {
            return T(2048) * t12 - T(6144) * t10 + T(6912) * t8 - T(3584) * t6 +
                   T(840) * t4 - T(72) * t2 + T(1);
          }
      }
    else
      {
        return cos(T(N) * atan2(ny, nx));
      }
  }

  /**
   * @brief Compute cos(N * psi) where theta = arctan(sqrt(nx^2+ny^2)/ny).
   *
   * This function provides an efficient way to calculate symmetries like $cos(N
   * arctan(\sqrt{n_x^2+n_y^2}/n_z))$, where `N` is a known number at compile time.
   *
   * Importantly, it doesn't always make sense to unroll the function with a
   * Chebyshev polynomial due to the crossover point in multiplication operations and
   * using the math library. We've found this to occur at N = 5, although your mileage
   * will vary based on your architecture.
   * TODO: Benchmark this better
   */
  template <unsigned int N, typename T>
  [[nodiscard]] inline DEAL_II_ALWAYS_INLINE T
  cos_psi(const T &nx, const T &ny, const T &nz) noexcept
  {
    using std::atan2;
    using std::cos;
    using std::sqrt;

    if constexpr (N <= 5)
      {
        const T t1  = nz;      // cos θ
        const T t2  = t1 * t1; // (cos θ)^2
        const T t3  = t2 * t1; // (cos θ)^3
        const T t4  = t2 * t2; // (cos θ)^4
        const T t5  = t3 * t2; // (cos θ)^5
        const T t6  = t3 * t3; // (cos θ)^6
        const T t7  = t4 * t3; // (cos θ)^7
        const T t8  = t4 * t4; // (cos θ)^8
        const T t9  = t5 * t4; // (cos θ)^9
        const T t10 = t5 * t5; // (cos θ)^10
        const T t11 = t6 * t5; // (cos θ)^11
        const T t12 = t6 * t6; // (cos θ)^12

        if constexpr (N == 0)
          {
            return T(1);
          }
        else if constexpr (N == 1)
          {
            return t1;
          }
        else if constexpr (N == 2)
          {
            return T(2) * t2 - T(1);
          }
        else if constexpr (N == 3)
          {
            return T(4) * t3 - T(3) * t1;
          }
        else if constexpr (N == 4)
          {
            return T(8) * t4 - T(8) * t2 + T(1);
          }
        else if constexpr (N == 5)
          {
            return T(16) * t5 - T(20) * t3 + T(5) * t1;
          }
        else if constexpr (N == 6)
          {
            return T(32) * t6 - T(48) * t4 + T(18) * t2 - T(1);
          }
        else if constexpr (N == 7)
          {
            return T(64) * t7 - T(112) * t5 + T(56) * t3 - T(7) * t1;
          }
        else if constexpr (N == 8)
          {
            return T(128) * t8 - T(256) * t6 + T(160) * t4 - T(32) * t2 + T(1);
          }
        else if constexpr (N == 9)
          {
            return T(256) * t9 - T(576) * t7 + T(432) * t5 - T(120) * t3 + T(9) * t1;
          }
        else if constexpr (N == 10)
          {
            return T(512) * t10 - T(1280) * t8 + T(1120) * t6 - T(400) * t4 + T(50) * t2 -
                   T(1);
          }
        else if constexpr (N == 11)
          {
            return T(1024) * t11 - T(2816) * t9 + T(2816) * t7 - T(1232) * t5 +
                   T(220) * t3 - T(11) * t1;
          }
        else if constexpr (N == 12)
          {
            return T(2048) * t12 - T(6144) * t10 + T(6912) * t8 - T(3584) * t6 +
                   T(840) * t4 - T(72) * t2 + T(1);
          }
      }
    else
      {
        return cos(T(N) * atan2(sqrt(nx * nx + ny * ny), nz));
      }
  }

} // namespace Symmetries

PRISMS_PF_END_NAMESPACE
