#include <deal.II/base/config.h>
#include <deal.II/base/vectorization.h>

#include <prismspf/utilities/crystal_symmetry.h>
#include <prismspf/utilities/vectorized_operations.h>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

template <unsigned int N, typename T>
inline DEAL_II_ALWAYS_INLINE T
reference_cos_arctan(const T &x)
{
  return std::cos(T(N) * std::atan(x));
}

template <unsigned int N, typename T>
inline DEAL_II_ALWAYS_INLINE T
reference_sin_arctan(const T &x)
{
  return std::sin(T(N) * std::atan(x));
}

template <unsigned int N, typename T>
inline DEAL_II_ALWAYS_INLINE T
reference_cos_theta(const T &nx, const T &ny)
{
  return std::cos(T(N) * std::atan2(ny, nx));
}

template <unsigned int N, typename T>
inline DEAL_II_ALWAYS_INLINE T
reference_sin_theta(const T &nx, const T &ny)
{
  return std::sin(T(N) * std::atan2(ny, nx));
}

template <unsigned int N, typename T>
inline DEAL_II_ALWAYS_INLINE T
reference_cos_psi(const T &nx, const T &ny, const T &nz)
{
  return std::cos(T(N) * std::atan2(std::sqrt(ny * ny + nx * nx), nz));
}

template <unsigned int N, typename T>
inline DEAL_II_ALWAYS_INLINE T
reference_sin_psi(const T &nx, const T &ny, const T &nz)
{
  return std::sin(T(N) * std::atan2(std::sqrt(ny * ny + nx * nx), nz));
}

template <typename Number, std::size_t width = dealii::VectorizedArray<Number>::size()>
void
check_all_lanes_abs(const dealii::VectorizedArray<Number, width> &result,
                    const dealii::VectorizedArray<Number, width> &expected,
                    const Number                                  tol)
{
  for (std::size_t i = 0; i < width; ++i)
    {
      CHECK_THAT(result[i], Catch::Matchers::WithinAbs(expected[i], tol));
    }
}

TEMPLATE_TEST_CASE("Symmetries::cos_arctan matches std::cos(N*std::atan(x))",
                   "[vectorized][symmetries][cos_arctan]",
                   float,
                   double)
{
  using VecArray     = dealii::VectorizedArray<TestType>;
  constexpr auto tol = TestType(1e-5);

  std::array<VecArray, 5> inputs {VecArray(0.0),
                                  VecArray(1.0),
                                  VecArray(-1.0),
                                  VecArray(10.0),
                                  VecArray(-10.0)};

  for (unsigned int N = 0; N <= 12; ++N)
    {
      DYNAMIC_SECTION("N = " << N)
      {
        for (const auto input : inputs)
          {
            VecArray value;
            VecArray expected;

            switch (N)
              {
                case 0:
                  value    = prisms::Symmetries::cos_arctan<0>(input);
                  expected = reference_cos_arctan<0>(input);
                  break;
                case 1:
                  value    = prisms::Symmetries::cos_arctan<1>(input);
                  expected = reference_cos_arctan<1>(input);
                  break;
                case 2:
                  value    = prisms::Symmetries::cos_arctan<2>(input);
                  expected = reference_cos_arctan<2>(input);
                  break;
                case 3:
                  value    = prisms::Symmetries::cos_arctan<3>(input);
                  expected = reference_cos_arctan<3>(input);
                  break;
                case 4:
                  value    = prisms::Symmetries::cos_arctan<4>(input);
                  expected = reference_cos_arctan<4>(input);
                  break;
                case 5:
                  value    = prisms::Symmetries::cos_arctan<5>(input);
                  expected = reference_cos_arctan<5>(input);
                  break;
                case 6:
                  value    = prisms::Symmetries::cos_arctan<6>(input);
                  expected = reference_cos_arctan<6>(input);
                  break;
                case 7:
                  value    = prisms::Symmetries::cos_arctan<7>(input);
                  expected = reference_cos_arctan<7>(input);
                  break;
                case 8:
                  value    = prisms::Symmetries::cos_arctan<8>(input);
                  expected = reference_cos_arctan<8>(input);
                  break;
                case 9:
                  value    = prisms::Symmetries::cos_arctan<9>(input);
                  expected = reference_cos_arctan<9>(input);
                  break;
                case 10:
                  value    = prisms::Symmetries::cos_arctan<10>(input);
                  expected = reference_cos_arctan<10>(input);
                  break;
                case 11:
                  value    = prisms::Symmetries::cos_arctan<11>(input);
                  expected = reference_cos_arctan<11>(input);
                  break;
                case 12:
                  value    = prisms::Symmetries::cos_arctan<12>(input);
                  expected = reference_cos_arctan<12>(input);
                  break;
              }

            check_all_lanes_abs(value, expected, tol);
          }
      }
    }
}

TEMPLATE_TEST_CASE("Symmetries::sin_arctan matches std::sin(N*std::atan(x))",
                   "[vectorized][symmetries][sin_arctan]",
                   float,
                   double)
{
  using VecArray     = dealii::VectorizedArray<TestType>;
  constexpr auto tol = TestType(1e-5);

  std::array<VecArray, 5> inputs {VecArray(0.0),
                                  VecArray(1.0),
                                  VecArray(-1.0),
                                  VecArray(10.0),
                                  VecArray(-10.0)};

  for (unsigned int N = 0; N <= 12; ++N)
    {
      DYNAMIC_SECTION("N = " << N)
      {
        for (const auto input : inputs)
          {
            VecArray value;
            VecArray expected;

            switch (N)
              {
                case 0:
                  value    = prisms::Symmetries::sin_arctan<0>(input);
                  expected = reference_sin_arctan<0>(input);
                  break;
                case 1:
                  value    = prisms::Symmetries::sin_arctan<1>(input);
                  expected = reference_sin_arctan<1>(input);
                  break;
                case 2:
                  value    = prisms::Symmetries::sin_arctan<2>(input);
                  expected = reference_sin_arctan<2>(input);
                  break;
                case 3:
                  value    = prisms::Symmetries::sin_arctan<3>(input);
                  expected = reference_sin_arctan<3>(input);
                  break;
                case 4:
                  value    = prisms::Symmetries::sin_arctan<4>(input);
                  expected = reference_sin_arctan<4>(input);
                  break;
                case 5:
                  value    = prisms::Symmetries::sin_arctan<5>(input);
                  expected = reference_sin_arctan<5>(input);
                  break;
                case 6:
                  value    = prisms::Symmetries::sin_arctan<6>(input);
                  expected = reference_sin_arctan<6>(input);
                  break;
                case 7:
                  value    = prisms::Symmetries::sin_arctan<7>(input);
                  expected = reference_sin_arctan<7>(input);
                  break;
                case 8:
                  value    = prisms::Symmetries::sin_arctan<8>(input);
                  expected = reference_sin_arctan<8>(input);
                  break;
                case 9:
                  value    = prisms::Symmetries::sin_arctan<9>(input);
                  expected = reference_sin_arctan<9>(input);
                  break;
                case 10:
                  value    = prisms::Symmetries::sin_arctan<10>(input);
                  expected = reference_sin_arctan<10>(input);
                  break;
                case 11:
                  value    = prisms::Symmetries::sin_arctan<11>(input);
                  expected = reference_sin_arctan<11>(input);
                  break;
                case 12:
                  value    = prisms::Symmetries::sin_arctan<12>(input);
                  expected = reference_sin_arctan<12>(input);
                  break;
              }

            check_all_lanes_abs(value, expected, tol);
          }
      }
    }
}

TEMPLATE_TEST_CASE("Symmetries::cos_theta matches std::cos(N*std::atan2(ny,nx))",
                   "[vectorized][symmetries][cos_theta]",
                   float,
                   double)
{
  using VecArray     = dealii::VectorizedArray<TestType>;
  constexpr auto tol = TestType(1e-5);

  std::array<std::pair<VecArray, VecArray>, 5> normals {
    {{VecArray(1.0), VecArray(0.0)},
     {VecArray(0.0), VecArray(1.0)},
     {VecArray(-1.0), VecArray(0.0)},
     {VecArray(1.0) / VecArray(std::sqrt(2.0)), VecArray(1.0) / VecArray(std::sqrt(2.0))},
     {VecArray(-1.0) / VecArray(std::sqrt(2.0)),
      VecArray(1.0) / VecArray(std::sqrt(2.0))}}
  };

  for (unsigned int N = 0; N <= 12; ++N)
    {
      DYNAMIC_SECTION("N = " << N)
      {
        for (const auto [nx, ny] : normals)
          {
            VecArray value;
            VecArray expected;

            switch (N)
              {
                case 0:
                  value    = prisms::Symmetries::cos_theta<0>(nx, ny);
                  expected = reference_cos_theta<0>(nx, ny);
                  break;
                case 1:
                  value    = prisms::Symmetries::cos_theta<1>(nx, ny);
                  expected = reference_cos_theta<1>(nx, ny);
                  break;
                case 2:
                  value    = prisms::Symmetries::cos_theta<2>(nx, ny);
                  expected = reference_cos_theta<2>(nx, ny);
                  break;
                case 3:
                  value    = prisms::Symmetries::cos_theta<3>(nx, ny);
                  expected = reference_cos_theta<3>(nx, ny);
                  break;
                case 4:
                  value    = prisms::Symmetries::cos_theta<4>(nx, ny);
                  expected = reference_cos_theta<4>(nx, ny);
                  break;
                case 5:
                  value    = prisms::Symmetries::cos_theta<5>(nx, ny);
                  expected = reference_cos_theta<5>(nx, ny);
                  break;
                case 6:
                  value    = prisms::Symmetries::cos_theta<6>(nx, ny);
                  expected = reference_cos_theta<6>(nx, ny);
                  break;
                case 7:
                  value    = prisms::Symmetries::cos_theta<7>(nx, ny);
                  expected = reference_cos_theta<7>(nx, ny);
                  break;
                case 8:
                  value    = prisms::Symmetries::cos_theta<8>(nx, ny);
                  expected = reference_cos_theta<8>(nx, ny);
                  break;
                case 9:
                  value    = prisms::Symmetries::cos_theta<9>(nx, ny);
                  expected = reference_cos_theta<9>(nx, ny);
                  break;
                case 10:
                  value    = prisms::Symmetries::cos_theta<10>(nx, ny);
                  expected = reference_cos_theta<10>(nx, ny);
                  break;
                case 11:
                  value    = prisms::Symmetries::cos_theta<11>(nx, ny);
                  expected = reference_cos_theta<11>(nx, ny);
                  break;
                case 12:
                  value    = prisms::Symmetries::cos_theta<12>(nx, ny);
                  expected = reference_cos_theta<12>(nx, ny);
                  break;
              }

            check_all_lanes_abs(value, expected, tol);
          }
      }
    }
}

TEMPLATE_TEST_CASE("Symmetries::sin_theta matches std::sin(N*std::atan2(ny,nx))",
                   "[vectorized][symmetries][sin_theta]",
                   float,
                   double)
{
  using VecArray     = dealii::VectorizedArray<TestType>;
  constexpr auto tol = TestType(1e-5);

  std::array<std::pair<VecArray, VecArray>, 5> normals {
    {{VecArray(1.0), VecArray(0.0)},
     {VecArray(0.0), VecArray(1.0)},
     {VecArray(-1.0), VecArray(0.0)},
     {VecArray(1.0) / VecArray(std::sqrt(2.0)), VecArray(1.0) / VecArray(std::sqrt(2.0))},
     {VecArray(-1.0) / VecArray(std::sqrt(2.0)),
      VecArray(1.0) / VecArray(std::sqrt(2.0))}}
  };

  for (unsigned int N = 0; N <= 12; ++N)
    {
      DYNAMIC_SECTION("N = " << N)
      {
        for (const auto [nx, ny] : normals)
          {
            VecArray value;
            VecArray expected;

            switch (N)
              {
                case 0:
                  value    = prisms::Symmetries::sin_theta<0>(nx, ny);
                  expected = reference_sin_theta<0>(nx, ny);
                  break;
                case 1:
                  value    = prisms::Symmetries::sin_theta<1>(nx, ny);
                  expected = reference_sin_theta<1>(nx, ny);
                  break;
                case 2:
                  value    = prisms::Symmetries::sin_theta<2>(nx, ny);
                  expected = reference_sin_theta<2>(nx, ny);
                  break;
                case 3:
                  value    = prisms::Symmetries::sin_theta<3>(nx, ny);
                  expected = reference_sin_theta<3>(nx, ny);
                  break;
                case 4:
                  value    = prisms::Symmetries::sin_theta<4>(nx, ny);
                  expected = reference_sin_theta<4>(nx, ny);
                  break;
                case 5:
                  value    = prisms::Symmetries::sin_theta<5>(nx, ny);
                  expected = reference_sin_theta<5>(nx, ny);
                  break;
                case 6:
                  value    = prisms::Symmetries::sin_theta<6>(nx, ny);
                  expected = reference_sin_theta<6>(nx, ny);
                  break;
                case 7:
                  value    = prisms::Symmetries::sin_theta<7>(nx, ny);
                  expected = reference_sin_theta<7>(nx, ny);
                  break;
                case 8:
                  value    = prisms::Symmetries::sin_theta<8>(nx, ny);
                  expected = reference_sin_theta<8>(nx, ny);
                  break;
                case 9:
                  value    = prisms::Symmetries::sin_theta<9>(nx, ny);
                  expected = reference_sin_theta<9>(nx, ny);
                  break;
                case 10:
                  value    = prisms::Symmetries::sin_theta<10>(nx, ny);
                  expected = reference_sin_theta<10>(nx, ny);
                  break;
                case 11:
                  value    = prisms::Symmetries::sin_theta<11>(nx, ny);
                  expected = reference_sin_theta<11>(nx, ny);
                  break;
                case 12:
                  value    = prisms::Symmetries::sin_theta<12>(nx, ny);
                  expected = reference_sin_theta<12>(nx, ny);
                  break;
              }

            check_all_lanes_abs(value, expected, tol);
          }
      }
    }
}

TEMPLATE_TEST_CASE(
  "Symmetries::cos_psi matches std::cos(N*std::atan2(sqrt(nx*nx+ny*ny),nz)",
  "[vectorized][symmetries][cos_psi]",
  float,
  double)
{
  using VecArray     = dealii::VectorizedArray<TestType>;
  constexpr auto tol = TestType(1e-5);

  std::array<std::tuple<VecArray, VecArray, VecArray>, 5> normals {
    {{VecArray(1.0), VecArray(0.0), VecArray(0.0)},
     {VecArray(0.0), VecArray(1.0), VecArray(0.0)},
     {VecArray(-1.0), VecArray(0.0), VecArray(0.0)},
     {VecArray(1.0) / VecArray(std::sqrt(3.0)),
      VecArray(1.0) / VecArray(std::sqrt(3.0)),
      VecArray(1.0) / VecArray(std::sqrt(3.0))},
     {VecArray(-1.0) / VecArray(std::sqrt(3.0)),
      VecArray(1.0) / VecArray(std::sqrt(3.0)),
      VecArray(1.0) / VecArray(std::sqrt(3.0))}}
  };

  for (unsigned int N = 0; N <= 12; ++N)
    {
      DYNAMIC_SECTION("N = " << N)
      {
        for (const auto &[nx, ny, nz] : normals)
          {
            VecArray value;
            VecArray expected;

            switch (N)
              {
                case 0:
                  value    = prisms::Symmetries::cos_psi<0>(nx, ny, nz);
                  expected = reference_cos_psi<0>(nx, ny, nz);
                  break;
                case 1:
                  value    = prisms::Symmetries::cos_psi<1>(nx, ny, nz);
                  expected = reference_cos_psi<1>(nx, ny, nz);
                  break;
                case 2:
                  value    = prisms::Symmetries::cos_psi<2>(nx, ny, nz);
                  expected = reference_cos_psi<2>(nx, ny, nz);
                  break;
                case 3:
                  value    = prisms::Symmetries::cos_psi<3>(nx, ny, nz);
                  expected = reference_cos_psi<3>(nx, ny, nz);
                  break;
                case 4:
                  value    = prisms::Symmetries::cos_psi<4>(nx, ny, nz);
                  expected = reference_cos_psi<4>(nx, ny, nz);
                  break;
                case 5:
                  value    = prisms::Symmetries::cos_psi<5>(nx, ny, nz);
                  expected = reference_cos_psi<5>(nx, ny, nz);
                  break;
                case 6:
                  value    = prisms::Symmetries::cos_psi<6>(nx, ny, nz);
                  expected = reference_cos_psi<6>(nx, ny, nz);
                  break;
                case 7:
                  value    = prisms::Symmetries::cos_psi<7>(nx, ny, nz);
                  expected = reference_cos_psi<7>(nx, ny, nz);
                  break;
                case 8:
                  value    = prisms::Symmetries::cos_psi<8>(nx, ny, nz);
                  expected = reference_cos_psi<8>(nx, ny, nz);
                  break;
                case 9:
                  value    = prisms::Symmetries::cos_psi<9>(nx, ny, nz);
                  expected = reference_cos_psi<9>(nx, ny, nz);
                  break;
                case 10:
                  value    = prisms::Symmetries::cos_psi<10>(nx, ny, nz);
                  expected = reference_cos_psi<10>(nx, ny, nz);
                  break;
                case 11:
                  value    = prisms::Symmetries::cos_psi<11>(nx, ny, nz);
                  expected = reference_cos_psi<11>(nx, ny, nz);
                  break;
                case 12:
                  value    = prisms::Symmetries::cos_psi<12>(nx, ny, nz);
                  expected = reference_cos_psi<12>(nx, ny, nz);
                  break;
              }

            check_all_lanes_abs(value, expected, tol);
          }
      }
    }
}

TEMPLATE_TEST_CASE(
  "Symmetries::sin_psi matches std::sin(N*std::atan2(sqrt(nx*nx+ny*ny),nz)",
  "[vectorized][symmetries][sin_psi]",
  float,
  double)
{
  using VecArray     = dealii::VectorizedArray<TestType>;
  constexpr auto tol = TestType(1e-5);

  std::array<std::tuple<VecArray, VecArray, VecArray>, 5> normals {
    {{VecArray(1.0), VecArray(0.0), VecArray(0.0)},
     {VecArray(0.0), VecArray(1.0), VecArray(0.0)},
     {VecArray(-1.0), VecArray(0.0), VecArray(0.0)},
     {VecArray(1.0) / VecArray(std::sqrt(3.0)),
      VecArray(1.0) / VecArray(std::sqrt(3.0)),
      VecArray(1.0) / VecArray(std::sqrt(3.0))},
     {VecArray(-1.0) / VecArray(std::sqrt(3.0)),
      VecArray(1.0) / VecArray(std::sqrt(3.0)),
      VecArray(1.0) / VecArray(std::sqrt(3.0))}}
  };

  for (unsigned int N = 0; N <= 12; ++N)
    {
      DYNAMIC_SECTION("N = " << N)
      {
        for (const auto &[nx, ny, nz] : normals)
          {
            VecArray value;
            VecArray expected;

            switch (N)
              {
                case 0:
                  value    = prisms::Symmetries::sin_psi<0>(nx, ny, nz);
                  expected = reference_sin_psi<0>(nx, ny, nz);
                  break;
                case 1:
                  value    = prisms::Symmetries::sin_psi<1>(nx, ny, nz);
                  expected = reference_sin_psi<1>(nx, ny, nz);
                  break;
                case 2:
                  value    = prisms::Symmetries::sin_psi<2>(nx, ny, nz);
                  expected = reference_sin_psi<2>(nx, ny, nz);
                  break;
                case 3:
                  value    = prisms::Symmetries::sin_psi<3>(nx, ny, nz);
                  expected = reference_sin_psi<3>(nx, ny, nz);
                  break;
                case 4:
                  value    = prisms::Symmetries::sin_psi<4>(nx, ny, nz);
                  expected = reference_sin_psi<4>(nx, ny, nz);
                  break;
                case 5:
                  value    = prisms::Symmetries::sin_psi<5>(nx, ny, nz);
                  expected = reference_sin_psi<5>(nx, ny, nz);
                  break;
                case 6:
                  value    = prisms::Symmetries::sin_psi<6>(nx, ny, nz);
                  expected = reference_sin_psi<6>(nx, ny, nz);
                  break;
                case 7:
                  value    = prisms::Symmetries::sin_psi<7>(nx, ny, nz);
                  expected = reference_sin_psi<7>(nx, ny, nz);
                  break;
                case 8:
                  value    = prisms::Symmetries::sin_psi<8>(nx, ny, nz);
                  expected = reference_sin_psi<8>(nx, ny, nz);
                  break;
                case 9:
                  value    = prisms::Symmetries::sin_psi<9>(nx, ny, nz);
                  expected = reference_sin_psi<9>(nx, ny, nz);
                  break;
                case 10:
                  value    = prisms::Symmetries::sin_psi<10>(nx, ny, nz);
                  expected = reference_sin_psi<10>(nx, ny, nz);
                  break;
                case 11:
                  value    = prisms::Symmetries::sin_psi<11>(nx, ny, nz);
                  expected = reference_sin_psi<11>(nx, ny, nz);
                  break;
                case 12:
                  value    = prisms::Symmetries::sin_psi<12>(nx, ny, nz);
                  expected = reference_sin_psi<12>(nx, ny, nz);
                  break;
              }

            check_all_lanes_abs(value, expected, tol);
          }
      }
    }
}
