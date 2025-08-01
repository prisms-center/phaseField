// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/config.h>
#include <prismspf/field_input/read_vtk.h>

#include "catch.hpp"

PRISMS_PF_BEGIN_NAMESPACE

TEST_CASE("Read vtk file")
{
  SECTION("2D Scalar fields with quadrature degree 1")
  {
    ReadUnstructuredVTK<2, double> reader("field_input/test_2D_1degree.vtk");

    // Check that the number of points and cells are correct
    REQUIRE(reader.get_n_points() == 16);
    REQUIRE(reader.get_n_cells() == 4);

    // Check that the points are correct
    for (double x_coord : {0.0, 5.0, 10.0})
      {
        for (double y_coord : {0.0, 5.0, 10.0})
          {
            dealii::Point<2> point = {x_coord, y_coord};
            REQUIRE(reader.get_scalar_value(point, "n") == x_coord);
            REQUIRE(reader.get_scalar_value(point, "|nx|") == 1.0);
          }
      }
  }
  SECTION("2D Scalar fields with quadrature degree 2")
  {
    ReadUnstructuredVTK<2, double> reader("field_input/test_2D_2degree.vtk");

    // Check that the number of points and cells are correct
    REQUIRE(reader.get_n_points() == 36);
    REQUIRE(reader.get_n_cells() == 16);

    // Check that the points are correct
    for (double x_coord : {0.0, 2.5, 5.0, 7.5, 10.0})
      {
        for (double y_coord : {0.0, 2.5, 5.0, 7.5, 10.0})
          {
            dealii::Point<2> point = {x_coord, y_coord};
            REQUIRE(reader.get_scalar_value(point, "n") == x_coord);
            REQUIRE(reader.get_scalar_value(point, "|nx|") == 1.0);
          }
      }
  }
  SECTION("3D Scalar fields with quadrature degree 1")
  {
    ReadUnstructuredVTK<3, double> reader("field_input/test_3D_1degree.vtk");

    // Check that the number of points and cells are correct
    REQUIRE(reader.get_n_points() == 64);
    REQUIRE(reader.get_n_cells() == 8);

    // Check that the points are correct
    for (double x_coord : {0.0, 5.0, 10.0})
      {
        for (double y_coord : {0.0, 5.0, 10.0})
          {
            for (double z_coord : {0.0, 5.0, 10.0})
              {
                dealii::Point<3> point = {x_coord, y_coord, z_coord};
                REQUIRE(reader.get_scalar_value(point, "n") == x_coord);
                REQUIRE(reader.get_scalar_value(point, "|nx|") == 1.0);
              }
          }
      }
  }
  SECTION("Nonexistent file")
  {
    REQUIRE_THROWS(ReadUnstructuredVTK<2, double>("nonexistent.vtk"));
  }
  SECTION("Invalid cell type due to mismatched dimensions")
  {
    REQUIRE_THROWS(ReadUnstructuredVTK<3, double>("field_input/test_2D_1degree.vtk"));
  }
}

PRISMS_PF_END_NAMESPACE