// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/tensor.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include "catch.hpp"

PRISMS_PF_BEGIN_NAMESPACE

/**
 * This unit test looks at the compute stress utility funcition
 */
TEST_CASE("Compute stress")
{
  SECTION("1D")
  {
    const unsigned int dim = 1;

    // Create objects
    dealii::Tensor<2, (2 * dim) - 1 + (dim / 3), double> elasticity_tensor;
    dealii::Tensor<1, (2 * dim) - 1 + (dim / 3), double> strain;
    dealii::Tensor<1, (2 * dim) - 1 + (dim / 3), double> stress;

    // Fill in elasticity tensor and strain
    elasticity_tensor[0][0] = 2.0;
    strain[0]               = 0.1;

    // Compute stress & validate
    compute_stress<dim, double>(elasticity_tensor, strain, stress);

    REQUIRE(elasticity_tensor.dimension == 1);
    REQUIRE(strain.dimension == 1);
    REQUIRE(stress.dimension == 1);
    REQUIRE(stress[0] == 0.2);
  }
  SECTION("2D")
  {
    const unsigned int dim = 2;

    // Create objects
    dealii::Tensor<2, (2 * dim) - 1 + (dim / 3), double> elasticity_tensor;
    dealii::Tensor<1, (2 * dim) - 1 + (dim / 3), double> strain;
    dealii::Tensor<1, (2 * dim) - 1 + (dim / 3), double> stress;

    // Fill in elasticity tensor and strain
    elasticity_tensor[0][0] = 2.0;
    elasticity_tensor[1][1] = 2.0;
    elasticity_tensor[0][1] = -0.1;
    elasticity_tensor[1][0] = -0.1;
    strain[0]               = 0.1;
    strain[1]               = 0.1;

    // Compute stress & validate
    compute_stress<dim, double>(elasticity_tensor, strain, stress);

    REQUIRE(elasticity_tensor.dimension == 3);
    REQUIRE(strain.dimension == 3);
    REQUIRE(stress.dimension == 3);
    REQUIRE(stress[0] == 0.19);
    REQUIRE(stress[1] == 0.19);
    REQUIRE(stress[2] == 0.0);
  }
  SECTION("3D")
  {
    const unsigned int dim = 3;

    // Create objects
    dealii::Tensor<2, (2 * dim) - 1 + (dim / 3), double> elasticity_tensor;
    dealii::Tensor<1, (2 * dim) - 1 + (dim / 3), double> strain;
    dealii::Tensor<1, (2 * dim) - 1 + (dim / 3), double> stress;

    // Fill in elasticity tensor and strain
    elasticity_tensor[0][0] = 2.0;
    elasticity_tensor[1][1] = 2.0;
    elasticity_tensor[2][2] = 2.0;
    strain[0]               = 0.1;
    strain[1]               = 0.1;
    strain[2]               = 0.1;

    // Compute stress & validate
    compute_stress<dim, double>(elasticity_tensor, strain, stress);

    REQUIRE(elasticity_tensor.dimension == 6);
    REQUIRE(strain.dimension == 6);
    REQUIRE(stress.dimension == 6);
    REQUIRE(stress[0] == 0.2);
    REQUIRE(stress[1] == 0.2);
    REQUIRE(stress[2] == 0.2);
  }
  SECTION("2D - C array")
  {
    const unsigned int dim       = 2;
    const double       tolerance = 1e-3;

    // Create objects
    dealii::Tensor<2, (2 * dim) - 1 + (dim / 3), double> elasticity_tensor;
    double                                               strain[dim][dim];
    double                                               stress[dim][dim];

    // Fill in elasticity tensor and strain
    elasticity_tensor[0][0] = 2.69231;
    elasticity_tensor[1][1] = 2.69231;
    elasticity_tensor[2][2] = 0.76923;
    elasticity_tensor[0][1] = 1.15385;
    elasticity_tensor[1][0] = 1.15385;
    strain[0][0]            = 1.0;
    strain[1][1]            = 0.0;
    strain[0][1]            = 0.0;
    strain[1][0]            = 0.0;

    // Compute stress & validate
    compute_stress<dim, double>(elasticity_tensor, strain, stress);

    REQUIRE((stress[0][0] - 2.69231) < tolerance);
    REQUIRE((stress[1][1] - 1.15385) < tolerance);
    REQUIRE((stress[0][1] - 0.0) < tolerance);
    REQUIRE((stress[1][0] - 0.0) < tolerance);
  }
}

PRISMS_PF_END_NAMESPACE
