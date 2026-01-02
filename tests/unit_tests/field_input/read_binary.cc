// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/types.h>

#include <prismspf/config.h>
#include <prismspf/field_input/read_binary.h>

#include "catch.hpp"

PRISMS_PF_BEGIN_NAMESPACE

TEST_CASE("Read binary file")
{
  // Number of points in each direction for all tests
  const unsigned int n_points_x = 11;
  const unsigned int n_points_y = 4;
  const unsigned int n_points_z = 3;

  // General InitialConditionFile structure that we modify for each test
  InitialConditionFile file {
    .filename                  = "",
    .dataset_format            = DataFormatType::FlatBinary,
    .file_variable_names       = {"n"},
    .simulation_variable_names = {"n"},
    .n_data_points             = {n_points_x, n_points_y, n_points_z}
  };

  const double mesh_size = 10.0;

  // Some step sizes
  const double dx = mesh_size / n_points_x;
  const double dy = mesh_size / n_points_y;
  const double dz = mesh_size / n_points_z;

  SECTION("1D Scalar field")
  {
    const unsigned int dim = 1;

    file.filename = "test_scalar_1D.bin";

    // Create a simple spatial discretization
    SpatialDiscretization<dim> spatial_discretization;
    for (unsigned int i : std::views::iota(0U, dim))
      {
        spatial_discretization.set_size(i, mesh_size);
      }
    spatial_discretization.postprocess_and_validate();

    // Grab the number of points and fill out a data array to write to file
    const unsigned int  total_points = n_points_x;
    std::vector<double> data(total_points);
    for (unsigned int i : std::views::iota(0U, total_points))
      {
        data[i] = i;
      }
    ReadBinary<dim, double>::write_file(data, file);

    // Read the file
    ReadBinary<dim, double> reader(file, spatial_discretization);

    // Check that the values are correct at each point and halfway between points
    for (unsigned int i : std::views::iota(0U, n_points_x))
      {
        const double       x_coord = i * dx;
        dealii::Point<dim> point(x_coord);

        // Evaluate at the grid point
        REQUIRE(reader.get_scalar_value(point, "n") ==
                Approx(i).epsilon(Defaults::tolerance));

        // Evaluate halfway between grid points (except at the last point)
        if (i < n_points_x - 1)
          {
            point[0] += 0.5 * dx;
            REQUIRE(reader.get_scalar_value(point, "n") ==
                    Approx(i + 0.5).epsilon(Defaults::tolerance));
          }
      }
  }
  SECTION("1D Vector field")
  {
    const unsigned int dim = 1;

    file.filename = "test_vector_1D.bin";

    // Create a simple spatial discretization
    SpatialDiscretization<dim> spatial_discretization;
    for (unsigned int i : std::views::iota(0U, dim))
      {
        spatial_discretization.set_size(i, mesh_size);
      }
    spatial_discretization.postprocess_and_validate();

    // Grab the number of points and fill out a data array to write to file
    const unsigned int  total_points = n_points_x;
    std::vector<double> data(dim * total_points);
    for (unsigned int i : std::views::iota(0U, total_points))
      {
        for (unsigned int d : std::views::iota(0U, dim))
          {
            data[(dim * i) + d] = i;
          }
      }
    ReadBinary<dim, double>::write_file(data, file);

    // Read the file
    ReadBinary<dim, double> reader(file, spatial_discretization);

    // Check that the values are correct at each point and halfway between points
    for (unsigned int i : std::views::iota(0U, n_points_x))
      {
        const double x_coord = i * dx;

        for (unsigned int d : std::views::iota(0U, dim))
          {
            dealii::Point<dim> point(x_coord);

            // Evaluate at the grid point
            REQUIRE(reader.get_vector_value(point, "n")[d] ==
                    Approx(i).epsilon(Defaults::tolerance));

            // Evaluate halfway between grid points (except at the last point)
            if (i < n_points_x - 1)
              {
                point[0] += 0.5 * dx;
                REQUIRE(reader.get_vector_value(point, "n")[d] ==
                        Approx(i + 0.5).epsilon(Defaults::tolerance));
              }
          }
      }
  }
  SECTION("2D Scalar field")
  {
    const unsigned int dim = 2;

    file.filename = "test_scalar_2D.bin";

    // Create a simple spatial discretization
    SpatialDiscretization<dim> spatial_discretization;
    for (unsigned int i : std::views::iota(0U, dim))
      {
        spatial_discretization.set_size(i, mesh_size);
      }
    spatial_discretization.postprocess_and_validate();

    // Grab the number of points and fill out a data array to write to file
    const unsigned int  total_points = n_points_x * n_points_y;
    std::vector<double> data(total_points);
    for (unsigned int i : std::views::iota(0U, total_points))
      {
        data[i] = i;
      }
    ReadBinary<dim, double>::write_file(data, file);

    // Read the file
    ReadBinary<dim, double> reader(file, spatial_discretization);

    // Check that the values are correct at each point and halfway between points
    for (unsigned int i : std::views::iota(0U, n_points_x))
      {
        const double x_coord = i * dx;

        for (unsigned int j : std::views::iota(0U, n_points_y))
          {
            const double y_coord = j * dy;

            dealii::Point<dim> point(x_coord, y_coord);

            // Evaluate at the grid point
            REQUIRE(reader.get_scalar_value(point, "n") ==
                    Approx(i + (n_points_x * j)).epsilon(Defaults::tolerance));

            // Evaluate halfway between grid points (except at the last point)
            if ((i < n_points_x - 1) && (j < n_points_y - 1))
              {
                // Halfway in x-direction
                point[0] += 0.5 * dx;
                REQUIRE(reader.get_scalar_value(point, "n") ==
                        Approx(i + 0.5 + (n_points_x * j)).epsilon(Defaults::tolerance));

                // Halfway in y-direction
                point[0] -= 0.5 * dx;
                point[1] += 0.5 * dy;
                REQUIRE(
                  reader.get_scalar_value(point, "n") ==
                  Approx(i + (n_points_x * (j + 0.5))).epsilon(Defaults::tolerance));

                // Halfway in both directions
                point[0] += 0.5 * dx;
                REQUIRE(reader.get_scalar_value(point, "n") ==
                        Approx(i + 0.5 + (n_points_x * (j + 0.5)))
                          .epsilon(Defaults::tolerance));
              }
          }
      }
  }
  SECTION("2D Vector field")
  {
    const unsigned int dim = 2;

    file.filename = "test_vector_2D.bin";

    // Create a simple spatial discretization
    SpatialDiscretization<dim> spatial_discretization;
    for (unsigned int i : std::views::iota(0U, dim))
      {
        spatial_discretization.set_size(i, mesh_size);
      }
    spatial_discretization.postprocess_and_validate();

    // Grab the number of points and fill out a data array to write to file
    const unsigned int  total_points = n_points_x * n_points_y;
    std::vector<double> data(dim * total_points);
    for (unsigned int i : std::views::iota(0U, total_points))
      {
        for (unsigned int d : std::views::iota(0U, dim))
          {
            data[(dim * i) + d] = i;
          }
      }
    ReadBinary<dim, double>::write_file(data, file);

    // Read the file
    ReadBinary<dim, double> reader(file, spatial_discretization);

    // Check that the values are correct at each point and halfway between points
    for (unsigned int i : std::views::iota(0U, n_points_x))
      {
        const double x_coord = i * dx;

        for (unsigned int j : std::views::iota(0U, n_points_y))
          {
            const double y_coord = j * dy;

            for (unsigned int d : std::views::iota(0U, dim))
              {
                dealii::Point<dim> point(x_coord, y_coord);

                // Evaluate at the grid point
                REQUIRE(reader.get_vector_value(point, "n")[d] ==
                        Approx(i + (n_points_x * j)).epsilon(Defaults::tolerance));

                // Evaluate halfway between grid points (except at the last point)
                if ((i < n_points_x - 1) && (j < n_points_y - 1))
                  {
                    // Halfway in x-direction
                    point[0] += 0.5 * dx;
                    REQUIRE(
                      reader.get_vector_value(point, "n")[d] ==
                      Approx(i + 0.5 + (n_points_x * j)).epsilon(Defaults::tolerance));

                    // Halfway in y-direction
                    point[0] -= 0.5 * dx;
                    point[1] += 0.5 * dy;
                    REQUIRE(
                      reader.get_vector_value(point, "n")[d] ==
                      Approx(i + (n_points_x * (j + 0.5))).epsilon(Defaults::tolerance));

                    // Halfway in both directions
                    point[0] += 0.5 * dx;
                    REQUIRE(reader.get_vector_value(point, "n")[d] ==
                            Approx(i + 0.5 + (n_points_x * (j + 0.5)))
                              .epsilon(Defaults::tolerance));
                  }
              }
          }
      }
  }
  SECTION("3D Scalar field")
  {
    const unsigned int dim = 3;

    file.filename = "test_scalar_3D.bin";

    // Create a simple spatial discretization
    SpatialDiscretization<dim> spatial_discretization;
    for (unsigned int i : std::views::iota(0U, dim))
      {
        spatial_discretization.set_size(i, mesh_size);
      }
    spatial_discretization.postprocess_and_validate();

    // Grab the number of points and fill out a data array to write to file
    const unsigned int  total_points = n_points_x * n_points_y * n_points_z;
    std::vector<double> data(total_points);
    for (unsigned int i : std::views::iota(0U, total_points))
      {
        data[i] = i;
      }
    ReadBinary<dim, double>::write_file(data, file);

    // Read the file
    ReadBinary<dim, double> reader(file, spatial_discretization);

    // Check that the values are correct at each point and halfway between points
    for (unsigned int i : std::views::iota(0U, n_points_x))
      {
        const double x_coord = i * dx;

        for (unsigned int j : std::views::iota(0U, n_points_y))
          {
            const double y_coord = j * dy;

            for (unsigned int k : std::views::iota(0U, n_points_z))
              {
                const double z_coord = k * dz;

                dealii::Point<dim> point(x_coord, y_coord, z_coord);

                // Evaluate at the grid point
                REQUIRE(reader.get_scalar_value(point, "n") ==
                        Approx(i + (n_points_x * j) + (n_points_x * n_points_y * k))
                          .epsilon(Defaults::tolerance));

                // Evaluate halfway between grid points (except at the last point)
                if ((i < n_points_x - 1) && (j < n_points_y - 1) && (k < n_points_z - 1))
                  {
                    // Halfway in x-direction
                    point[0] += 0.5 * dx;
                    REQUIRE(
                      reader.get_scalar_value(point, "n") ==
                      Approx(i + 0.5 + (n_points_x * j) + (n_points_x * n_points_y * k))
                        .epsilon(Defaults::tolerance));

                    // Halfway in y-direction
                    point[0] -= 0.5 * dx;
                    point[1] += 0.5 * dy;
                    REQUIRE(
                      reader.get_scalar_value(point, "n") ==
                      Approx(i + (n_points_x * (j + 0.5)) + (n_points_x * n_points_y * k))
                        .epsilon(Defaults::tolerance));

                    // Halfway in z-direction
                    point[1] -= 0.5 * dy;
                    point[2] += 0.5 * dz;
                    REQUIRE(
                      reader.get_scalar_value(point, "n") ==
                      Approx(i + (n_points_x * j) + (n_points_x * n_points_y * (k + 0.5)))
                        .epsilon(Defaults::tolerance));

                    // Halfway in all directions
                    point[0] += 0.5 * dx;
                    point[1] += 0.5 * dy;
                    REQUIRE(reader.get_scalar_value(point, "n") ==
                            Approx(i + 0.5 + (n_points_x * (j + 0.5)) +
                                   (n_points_x * n_points_y * (k + 0.5)))
                              .epsilon(Defaults::tolerance));
                  }
              }
          }
      }
  }
  SECTION("3D Vector field")
  {
    const unsigned int dim = 3;

    file.filename = "test_vector_3D.bin";

    // Create a simple spatial discretization
    SpatialDiscretization<dim> spatial_discretization;
    for (unsigned int i : std::views::iota(0U, dim))
      {
        spatial_discretization.set_size(i, mesh_size);
      }
    spatial_discretization.postprocess_and_validate();

    // Grab the number of points and fill out a data array to write to file
    const unsigned int  total_points = n_points_x * n_points_y * n_points_z;
    std::vector<double> data(dim * total_points);
    for (unsigned int i : std::views::iota(0U, total_points))
      {
        for (unsigned int d : std::views::iota(0U, dim))
          {
            data[(dim * i) + d] = i;
          }
      }
    ReadBinary<dim, double>::write_file(data, file);

    // Read the file
    ReadBinary<dim, double> reader(file, spatial_discretization);

    // Check that the values are correct at each point and halfway between points
    for (unsigned int i : std::views::iota(0U, n_points_x))
      {
        const double x_coord = i * dx;

        for (unsigned int j : std::views::iota(0U, n_points_y))
          {
            const double y_coord = j * dy;

            for (unsigned int k : std::views::iota(0U, n_points_z))
              {
                const double z_coord = k * dz;

                for (unsigned int d : std::views::iota(0U, dim))
                  {
                    dealii::Point<dim> point(x_coord, y_coord, z_coord);

                    // Evaluate at the grid point
                    REQUIRE(reader.get_vector_value(point, "n")[d] ==
                            Approx(i + (n_points_x * j) + (n_points_x * n_points_y * k))
                              .epsilon(Defaults::tolerance));

                    // Evaluate halfway between grid points (except at the last point)
                    if ((i < n_points_x - 1) && (j < n_points_y - 1) &&
                        (k < n_points_z - 1))
                      {
                        // Halfway in x-direction
                        point[0] += 0.5 * dx;
                        REQUIRE(reader.get_vector_value(point, "n")[d] ==
                                Approx(i + 0.5 + (n_points_x * j) +
                                       (n_points_x * n_points_y * k))
                                  .epsilon(Defaults::tolerance));

                        // Halfway in y-direction
                        point[0] -= 0.5 * dx;
                        point[1] += 0.5 * dy;
                        REQUIRE(reader.get_vector_value(point, "n")[d] ==
                                Approx(i + (n_points_x * (j + 0.5)) +
                                       (n_points_x * n_points_y * k))
                                  .epsilon(Defaults::tolerance));

                        // Halfway in z-direction
                        point[1] -= 0.5 * dy;
                        point[2] += 0.5 * dz;
                        REQUIRE(reader.get_vector_value(point, "n")[d] ==
                                Approx(i + (n_points_x * j) +
                                       (n_points_x * n_points_y * (k + 0.5)))
                                  .epsilon(Defaults::tolerance));

                        // Halfway in all directions
                        point[0] += 0.5 * dx;
                        point[1] += 0.5 * dy;
                        REQUIRE(reader.get_vector_value(point, "n")[d] ==
                                Approx(i + 0.5 + (n_points_x * (j + 0.5)) +
                                       (n_points_x * n_points_y * (k + 0.5)))
                                  .epsilon(Defaults::tolerance));
                      }
                  }
              }
          }
      }
  }
  SECTION("Nonexistent file")
  {
    file.filename = "nonexistent.bin";

    // Create a simple spatial discretization
    SpatialDiscretization<1> spatial_discretization;
    spatial_discretization.set_size(0, mesh_size);
    spatial_discretization.postprocess_and_validate();

    REQUIRE_THROWS(ReadBinary<1, double>(file, spatial_discretization));
  }
  SECTION("Invalid mesh type")
  {
    file.filename = "invalid_mesh_type.bin";

    // Create a simple spatial discretization
    SpatialDiscretization<1> spatial_discretization;
    spatial_discretization.set_radius(mesh_size);
    spatial_discretization.postprocess_and_validate();

    REQUIRE_THROWS(ReadBinary<1, double>(file, spatial_discretization));
  }
}

PRISMS_PF_END_NAMESPACE
