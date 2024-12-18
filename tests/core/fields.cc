#include "catch.hpp"

#include <core/fields.h>

/**
 * This unit test looks at fields.h and the initialization of various field types (PDEtype
 * & dim & spacedim). For successive objects of the same template parameters, the field
 * index should increase by 1. An error should be thrown for invalid enum types.
 */
TEST_CASE("Field declarations")
{
  const std::string field_name = "phi";

  SECTION("SCALAR")
  {
    SECTION("EXPLICIT_TIME_DEPENDENT")
    {
      Field<1> field_1d(SCALAR, EXPLICIT_TIME_DEPENDENT, field_name);
      Field<2> field_2d(SCALAR, EXPLICIT_TIME_DEPENDENT, field_name);
      Field<3> field_3d(SCALAR, EXPLICIT_TIME_DEPENDENT, field_name);

      REQUIRE(field_1d.numComponents == 1);
      REQUIRE(field_1d.index == 0);

      REQUIRE(field_2d.numComponents == 1);
      REQUIRE(field_2d.index == 0);

      REQUIRE(field_3d.numComponents == 1);
      REQUIRE(field_3d.index == 0);
    }
    SECTION("IMPLICIT_TIME_DEPENDENT")
    {
      Field<1> field_1d(SCALAR, IMPLICIT_TIME_DEPENDENT, field_name);
      Field<2> field_2d(SCALAR, IMPLICIT_TIME_DEPENDENT, field_name);
      Field<3> field_3d(SCALAR, IMPLICIT_TIME_DEPENDENT, field_name);

      REQUIRE(field_1d.numComponents == 1);
      REQUIRE(field_1d.index == 1);

      REQUIRE(field_2d.numComponents == 1);
      REQUIRE(field_2d.index == 1);

      REQUIRE(field_3d.numComponents == 1);
      REQUIRE(field_3d.index == 1);
    }
    SECTION("TIME_INDEPENDENT")
    {
      Field<1> field_1d(SCALAR, TIME_INDEPENDENT, field_name);
      Field<2> field_2d(SCALAR, TIME_INDEPENDENT, field_name);
      Field<3> field_3d(SCALAR, TIME_INDEPENDENT, field_name);

      REQUIRE(field_1d.numComponents == 1);
      REQUIRE(field_1d.index == 2);

      REQUIRE(field_2d.numComponents == 1);
      REQUIRE(field_2d.index == 2);

      REQUIRE(field_3d.numComponents == 1);
      REQUIRE(field_3d.index == 2);
    }
    SECTION("AUXILIARY")
    {
      Field<1> field_1d(SCALAR, AUXILIARY, field_name);
      Field<2> field_2d(SCALAR, AUXILIARY, field_name);
      Field<3> field_3d(SCALAR, AUXILIARY, field_name);

      REQUIRE(field_1d.numComponents == 1);
      REQUIRE(field_1d.index == 3);

      REQUIRE(field_2d.numComponents == 1);
      REQUIRE(field_2d.index == 3);

      REQUIRE(field_3d.numComponents == 1);
      REQUIRE(field_3d.index == 3);
    }
  }
  SECTION("VECTOR")
  {
    SECTION("EXPLICIT_TIME_DEPENDENT")
    {
      Field<1> field_1d(VECTOR, EXPLICIT_TIME_DEPENDENT, field_name);
      Field<2> field_2d(VECTOR, EXPLICIT_TIME_DEPENDENT, field_name);
      Field<3> field_3d(VECTOR, EXPLICIT_TIME_DEPENDENT, field_name);

      REQUIRE(field_1d.numComponents == 1);
      REQUIRE(field_1d.index == 4);

      REQUIRE(field_2d.numComponents == 2);
      REQUIRE(field_2d.index == 4);

      REQUIRE(field_3d.numComponents == 3);
      REQUIRE(field_3d.index == 4);
    }
    SECTION("IMPLICIT_TIME_DEPENDENT")
    {
      Field<1> field_1d(VECTOR, IMPLICIT_TIME_DEPENDENT, field_name);
      Field<2> field_2d(VECTOR, IMPLICIT_TIME_DEPENDENT, field_name);
      Field<3> field_3d(VECTOR, IMPLICIT_TIME_DEPENDENT, field_name);

      REQUIRE(field_1d.numComponents == 1);
      REQUIRE(field_1d.index == 5);

      REQUIRE(field_2d.numComponents == 2);
      REQUIRE(field_2d.index == 5);

      REQUIRE(field_3d.numComponents == 3);
      REQUIRE(field_3d.index == 5);
    }
    SECTION("TIME_INDEPENDENT")
    {
      Field<1> field_1d(VECTOR, TIME_INDEPENDENT, field_name);
      Field<2> field_2d(VECTOR, TIME_INDEPENDENT, field_name);
      Field<3> field_3d(VECTOR, TIME_INDEPENDENT, field_name);

      REQUIRE(field_1d.numComponents == 1);
      REQUIRE(field_1d.index == 6);

      REQUIRE(field_2d.numComponents == 2);
      REQUIRE(field_2d.index == 6);

      REQUIRE(field_3d.numComponents == 3);
      REQUIRE(field_3d.index == 6);
    }
    SECTION("AUXILIARY")
    {
      Field<1> field_1d(VECTOR, AUXILIARY, field_name);
      Field<2> field_2d(VECTOR, AUXILIARY, field_name);
      Field<3> field_3d(VECTOR, AUXILIARY, field_name);

      REQUIRE(field_1d.numComponents == 1);
      REQUIRE(field_1d.index == 7);

      REQUIRE(field_2d.numComponents == 2);
      REQUIRE(field_2d.index == 7);

      REQUIRE(field_3d.numComponents == 3);
      REQUIRE(field_3d.index == 7);
    }
  }
  SECTION("Invalid field type")
  {
    const auto invalid_field_type = static_cast<fieldType>(-1);

    REQUIRE_THROWS_AS(Field<3>(invalid_field_type, EXPLICIT_TIME_DEPENDENT, field_name),
                      std::invalid_argument);
  }

  SECTION("Invalid PDE type")
  {
    const auto invalid_PDE_type = static_cast<PDEType>(-1);

    REQUIRE_THROWS_AS(Field<3>(SCALAR, invalid_PDE_type, field_name),
                      std::invalid_argument);
  }
}