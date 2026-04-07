#include <deal.II/base/vectorization.h>

#include <prismspf/utilities/vectorized_operations.h>

#include <catch2/catch_template_test_macros.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>

template <typename Number, std::size_t width = dealii::VectorizedArray<Number>::size()>
void
check_all_lanes_rel(const dealii::VectorizedArray<Number, width> &result,
                    const Number                                  expected,
                    const Number                                  tol)
{
  for (std::size_t i = 0; i < width; ++i)
    {
      CHECK_THAT(result[i], Catch::Matchers::WithinRel(expected, tol));
    }
}

TEMPLATE_TEST_CASE("erf matches std::erf element-wise",
                   "[vectorized][erf]",
                   float,
                   double)
{
  using VecArray              = dealii::VectorizedArray<TestType>;
  constexpr std::size_t width = VecArray::size();

  constexpr auto tol = TestType(1e-6);

  SECTION("known values")
  {
    const std::array<TestType, 4> inputs = {TestType(0.0),
                                            TestType(1.0),
                                            TestType(-1.0),
                                            TestType(0.5)};

    for (const auto val : inputs)
      {
        VecArray x = val;
        check_all_lanes_rel(std::erf(x), std::erf(val), tol);
      }
  }
}

TEMPLATE_TEST_CASE("erfc matches std::erfc element-wise",
                   "[vectorized][erfc]",
                   float,
                   double)
{
  using VecArray              = dealii::VectorizedArray<TestType>;
  constexpr std::size_t width = VecArray::size();

  constexpr auto tol = TestType(1e-6);

  SECTION("known values")
  {
    const std::array<TestType, 4> inputs = {TestType(0.0),
                                            TestType(1.0),
                                            TestType(-1.0),
                                            TestType(0.5)};

    for (const auto val : inputs)
      {
        VecArray x = val;
        check_all_lanes_rel(std::erfc(x), std::erfc(val), tol);
      }
  }
}

TEMPLATE_TEST_CASE("fmod matches std::fmod element-wise",
                   "[vectorized][fmod]",
                   float,
                   double)
{
  using VecArray              = dealii::VectorizedArray<TestType>;
  constexpr std::size_t width = VecArray::size();

  constexpr auto tol = TestType(1e-6);

  SECTION("known values - single denom")
  {
    const std::array<std::pair<TestType, TestType>, 4> inputs = {
      std::pair {TestType(5.0),  TestType(3.0)},
      std::pair {TestType(1.0),  TestType(1.0)},
      std::pair {TestType(-5.0), TestType(3.0)},
      std::pair {TestType(0.0),  TestType(3.0)},
    };

    for (const auto &[dividend, divisor] : inputs)
      {
        VecArray x = dividend;
        TestType y = divisor;
        check_all_lanes_rel(std::fmod(x, y), std::fmod(dividend, divisor), tol);
      }
  }

  SECTION("known values - vectorized denom")
  {
    const std::array<std::pair<TestType, TestType>, 4> inputs = {
      std::pair {TestType(5.0),  TestType(3.0)},
      std::pair {TestType(1.0),  TestType(1.0)},
      std::pair {TestType(-5.0), TestType(3.0)},
      std::pair {TestType(0.0),  TestType(3.0)},
    };

    for (const auto &[dividend, divisor] : inputs)
      {
        VecArray x = dividend;
        VecArray y = divisor;
        check_all_lanes_rel(std::fmod(x, y), std::fmod(dividend, divisor), tol);
      }
  }
}
