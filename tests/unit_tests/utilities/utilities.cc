// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/tensor.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#include "catch.hpp"

PRISMS_PF_BEGIN_NAMESPACE

/**
 * This unit test looks at the string to type utility functions
 */
TEST_CASE("String to types")
{
  const std::unordered_map<std::string, unsigned int> unsigned_int_map = {
    {"",    0},
    {"one", 1},
    {"two", 2}
  };
  const std::unordered_map<std::string, int> int_map = {
    {"one",      1 },
    {"negative", -1}
  };

  SECTION("Single")
  {
    REQUIRE(string_to_type("one", unsigned_int_map) == 1);
    REQUIRE_THROWS(string_to_type("three", int_map));
  }
  SECTION("Pair")
  {
    auto pair_one = string_to_type_pair_with_delimiters("negative",
                                                        unsigned_int_map,
                                                        int_map,
                                                        {'(', ')'});
    auto pair_two = string_to_type_pair_with_delimiters("two(negative)",
                                                        unsigned_int_map,
                                                        int_map,
                                                        {'(', ')'});

    REQUIRE(pair_one.first == -1);
    REQUIRE(pair_one.second == 0);
    REQUIRE(pair_two.first == -1);
    REQUIRE(pair_two.second == 2);

    REQUIRE_THROWS(string_to_type_pair_with_delimiters("invalid",
                                                       unsigned_int_map,
                                                       int_map,
                                                       {'(', ')'}));
    REQUIRE_THROWS(string_to_type_pair_with_delimiters("invalid(negative)",
                                                       unsigned_int_map,
                                                       int_map,
                                                       {'(', ')'}));
    REQUIRE_THROWS(string_to_type_pair_with_delimiters("two(negative)(extra)",
                                                       unsigned_int_map,
                                                       int_map,
                                                       {'(', ')'}));
    REQUIRE_THROWS(string_to_type_pair_with_delimiters("two)(negative)",
                                                       unsigned_int_map,
                                                       int_map,
                                                       {'(', ')'}));
    REQUIRE_THROWS(string_to_type_pair_with_delimiters("two((negative))",
                                                       unsigned_int_map,
                                                       int_map,
                                                       {'(', ')'}));
    REQUIRE_THROWS(string_to_type_pair_with_delimiters("two{negative}",
                                                       unsigned_int_map,
                                                       int_map,
                                                       {'(', ')'}));
  }
}

/**
 * This unit test looks at the compute stress utility function
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
}

PRISMS_PF_END_NAMESPACE
