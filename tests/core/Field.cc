#include "../../contrib/catch/catch.hpp"
#include "../../include/fields.h"

TEST_CASE("Field declaration in 1D")
{
  const std::string field_name = "phi";

  SECTION("SCALAR")
  {
    SECTION("EXPLICIT_TIME_DEPENDENT")
    {
      Field<1> field(SCALAR, EXPLICIT_TIME_DEPENDENT, field_name);

      REQUIRE(field.numComponents == 1);
      REQUIRE(field.index == 0);
    }
    SECTION("IMPLICIT_TIME_DEPENDENT")
    {
      Field<1> field(SCALAR, IMPLICIT_TIME_DEPENDENT, field_name);

      REQUIRE(field.numComponents == 1);
      REQUIRE(field.index == 1);
    }
    SECTION("TIME_INDEPENDENT")
    {
      Field<1> field(SCALAR, TIME_INDEPENDENT, field_name);

      REQUIRE(field.numComponents == 1);
      REQUIRE(field.index == 2);
    }
    SECTION("AUXILIARY")
    {
      Field<1> field(SCALAR, AUXILIARY, field_name);

      REQUIRE(field.numComponents == 1);
      REQUIRE(field.index == 3);
    }
  }
  SECTION("VECTOR")
  {
    SECTION("EXPLICIT_TIME_DEPENDENT")
    {
      Field<1> field(SCALAR, EXPLICIT_TIME_DEPENDENT, field_name);

      REQUIRE(field.numComponents == 1);
      REQUIRE(field.index == 4);
    }
    SECTION("IMPLICIT_TIME_DEPENDENT")
    {
      Field<1> field(SCALAR, IMPLICIT_TIME_DEPENDENT, field_name);

      REQUIRE(field.numComponents == 1);
      REQUIRE(field.index == 5);
    }
    SECTION("TIME_INDEPENDENT")
    {
      Field<1> field(SCALAR, TIME_INDEPENDENT, field_name);

      REQUIRE(field.numComponents == 1);
      REQUIRE(field.index == 6);
    }
    SECTION("AUXILIARY")
    {
      Field<1> field(SCALAR, AUXILIARY, field_name);

      REQUIRE(field.numComponents == 1);
      REQUIRE(field.index == 7);
    }
  }
}