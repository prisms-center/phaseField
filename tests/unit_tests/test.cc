#include "../../src/models/mechanics/computeStress.h"

#define CATCH_CONFIG_MAIN
#include "../../contrib/catch/catch.hpp"

TEST_CASE("Compute Stress")
{
  SECTION("1D")
  {
    REQUIRE(4 == 4);
  }
  SECTION("2D")
  {
    REQUIRE(1 == 4);
  }
  SECTION("3D")
  {
    REQUIRE(1 == 4);
  }
}