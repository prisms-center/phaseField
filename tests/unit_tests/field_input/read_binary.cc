// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/config.h>
#include <prismspf/field_input/read_binary.h>

#include "catch.hpp"

PRISMS_PF_BEGIN_NAMESPACE



TEST_CASE("Read binary file")
{
  SECTION("Nonexistent file")
  {
    InitialConditionFile file;
    file.filename = "nonexistent.bin";
    SpatialDiscretization<2> spatial_discretization;
    REQUIRE_THROWS(ReadBinary<2, double>(file, spatial_discretization));
  }
  SECTION("1D read/write")
  {
    InitialConditionFile file;
    std::vector<double> data = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
    file.filename = "field_input/test_1D.bin";
    file.dataset_format = DataFormatType::FlatBinary;
    file.n_data_points = {6, 0, 0};
    ReadBinary<1, double>::write_file(data, file);
    SpatialDiscretization<1> spatial_discretization;
    spatial_discretization.set_size(0, 5.0);
    ReadBinary<1, double> reader(file, spatial_discretization);
    // Check values at selected points
    dealii::Point<1> point1(3.0);
    REQUIRE(reader.get_scalar_value(point1, "n") - 3.0 < 1e-12);
    dealii::Point<1> point2(2.5);
    REQUIRE(reader.get_scalar_value(point2, "n") - 2.5 < 1e-12);
  }
  SECTION("2D read/write")
  {
    InitialConditionFile file;
    std::vector<double> data = {0.0, 0.0, 0.0, 0.0,
                                0.0, 10.0, 0.0, 1.0,
                                0.0, 0.0, 0.0, 0.0};
    file.filename = "field_input/test_2D.bin";
    file.dataset_format = DataFormatType::FlatBinary;
    file.n_data_points = {4, 3, 0};
    ReadBinary<2, double>::write_file(data, file);
    SpatialDiscretization<2> spatial_discretization;
    spatial_discretization.set_size(0, 3.0);
    spatial_discretization.set_size(1, 2.0);
    ReadBinary<2, double> reader(file, spatial_discretization);

    // Check values at selected points
    dealii::Point<2> point1 = {0.5, 0.5};
    REQUIRE(reader.get_scalar_value(point1, "n") - 2.5 < 1e-12);
    dealii::Point<2> point2 = {3.0, 1.5};
    REQUIRE(reader.get_scalar_value(point2, "n") - 0.5 < 1e-12);
  }
  SECTION("3D read/write")
  {
    InitialConditionFile file;
    std::vector<double> data = {0.0, 0.0, 1.0, 1.0,
                                2.0, 2.0, 3.0, 3.0};
    file.filename = "field_input/test_3D.bin";
    file.dataset_format = DataFormatType::FlatBinary;
    file.n_data_points = {2, 2, 2};
    ReadBinary<3, double>::write_file(data, file);
    SpatialDiscretization<3> spatial_discretization;
    spatial_discretization.set_size(0, 1.0);
    spatial_discretization.set_size(1, 1.0);
    spatial_discretization.set_size(2, 1.0);
    ReadBinary<3, double> reader(file, spatial_discretization);

    // Check values at selected points
    dealii::Point<3> point1 = {0.5, 0.5, 0.2};
    REQUIRE(reader.get_scalar_value(point1, "n") - 0.9 < 1e-12);
    dealii::Point<3> point2 = {1.0, 1.0, 0.0};
    REQUIRE(reader.get_scalar_value(point2, "n") - 3.0 < 1e-12);
    dealii::Point<3> point3 = {0.0, 0.5, 1.0};
    REQUIRE(reader.get_scalar_value(point3, "n") - 2.5 < 1e-12);
  }
  SECTION("3D vector read/write")
  {
    InitialConditionFile file;
    std::vector<double> data = {0.0, 10.0, 20.0, 0.0, 10.0, 20.0,
                                1.0, 11.0, 21.0, 1.0, 11.0, 21.0,
                                2.0, 12.0, 22.0, 2.0, 12.0, 22.0,
                                3.0, 13.0, 23.0, 3.0, 13.0, 23.0};
    file.filename = "field_input/test_3D_vec.bin";
    file.dataset_format = DataFormatType::FlatBinary;
    file.n_data_points = {2, 2, 2};
    ReadBinary<3, double>::write_file(data, file);
    SpatialDiscretization<3> spatial_discretization;
    spatial_discretization.set_size(0, 1.0);
    spatial_discretization.set_size(1, 1.0);
    spatial_discretization.set_size(2, 1.0);
    ReadBinary<3, double> reader(file, spatial_discretization);

    // Check values at selected points
    dealii::Point<3> point1 = {0.5, 0.5, 0.2};
    REQUIRE(reader.get_vector_value(point1, "n")[2] - 20.9 < 1e-12);
    dealii::Point<3> point2 = {0.0, 0.0, 0.0};
    REQUIRE(reader.get_vector_value(point2, "n")[1] - 10.0 < 1e-12);
    dealii::Point<3> point3 = {0.0, 0.5, 1.0};
    REQUIRE(reader.get_vector_value(point3, "n")[0] - 2.5 < 1e-12);
  }
  SECTION("trigger warning for size mismatch")
  {
    InitialConditionFile file;
    file.filename = "field_input/test_3D_vec.bin"; // from previous test
    file.dataset_format = DataFormatType::FlatBinary;
    file.n_data_points = {20, 2, 0};
    SpatialDiscretization<3> spatial_discretization;
    REQUIRE_THROWS(ReadBinary<3, double>(file, spatial_discretization));
  }
}

PRISMS_PF_END_NAMESPACE