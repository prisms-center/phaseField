// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/utilities/mechanics.h>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <initializer_list>
#include <numbers>
#include <utility>

namespace
{
  namespace Mechanics = prismspf::Mechanics;

  using prismspf::StressState;

  constexpr StressState three_dimensional = StressState::ThreeDimensional;
  constexpr StressState plane_stress      = StressState::PlaneStress;
  constexpr StressState plane_strain      = StressState::PlaneStrain;

  constexpr double tolerance = 1.0e-12;

  void
  check_close(const double actual, const double expected)
  {
    CHECK_THAT(actual, Catch::Matchers::WithinAbs(expected, tolerance));
  }

  template <unsigned int size>
  dealii::Tensor<1, size, double>
  make_voigt(const std::initializer_list<double> values)
  {
    REQUIRE(values.size() == size);

    dealii::Tensor<1, size, double> result;

    unsigned int index = 0;
    for (const double value : values)
      result[index++] = value;

    return result;
  }

  template <unsigned int size>
  void
  check_voigt(const dealii::Tensor<1, size, double> &actual,
              const std::initializer_list<double>    expected)
  {
    REQUIRE(expected.size() == size);

    unsigned int index = 0;
    for (const double value : expected)
      check_close(actual[index++], value);
  }

  template <unsigned int dim, typename ActualTensor, typename ExpectedTensor>
  void
  check_rank_2_tensor(const ActualTensor &actual, const ExpectedTensor &expected)
  {
    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = 0; j < dim; ++j)
        check_close(actual[i][j], expected[i][j]);
  }

  template <unsigned int size>
  void
  check_matrix(const dealii::Tensor<2, size, double> &actual,
               const dealii::Tensor<2, size, double> &expected)
  {
    check_rank_2_tensor<size>(actual, expected);
  }

} // namespace

TEST_CASE("Mechanics voigt_size", "[mechanics][voigt][size]")
{
  STATIC_REQUIRE(Mechanics::voigt_size<1> == 1);
  STATIC_REQUIRE(Mechanics::voigt_size<3> == 6);

  STATIC_REQUIRE((Mechanics::voigt_size<2, StressState::PlaneStress> == 3));
  STATIC_REQUIRE((Mechanics::voigt_size<2, StressState::PlaneStrain> == 4));

  STATIC_REQUIRE((Mechanics::valid_stress_state<1, StressState::ThreeDimensional>) );
  STATIC_REQUIRE((Mechanics::valid_stress_state<2, StressState::PlaneStress>) );
  STATIC_REQUIRE((Mechanics::valid_stress_state<2, StressState::PlaneStrain>) );
  STATIC_REQUIRE((Mechanics::valid_stress_state<3, StressState::ThreeDimensional>) );

  STATIC_REQUIRE_FALSE(
    (Mechanics::valid_stress_state<2, StressState::ThreeDimensional>) );
  STATIC_REQUIRE_FALSE((Mechanics::valid_stress_state<3, StressState::PlaneStress>) );
}

TEST_CASE("Mechanics strain_to_voigt", "[mechanics][strain][voigt]")
{
  SECTION("1D Tensor")
  {
    dealii::Tensor<2, 1, double> strain;
    strain[0][0] = 2.5;

    const auto returned = Mechanics::strain_to_voigt<1>(strain);

    dealii::Tensor<1, 1, double> output;
    Mechanics::strain_to_voigt<1>(strain, output);

    check_voigt<1>(returned, {2.5});
    check_voigt<1>(output, {2.5});
  }

  SECTION("1D SymmetricTensor")
  {
    dealii::SymmetricTensor<2, 1, double> strain;
    strain[0][0] = 2.5;

    const auto returned = Mechanics::strain_to_voigt<1>(strain);

    dealii::Tensor<1, 1, double> output;
    Mechanics::strain_to_voigt<1>(strain, output);

    check_voigt<1>(returned, {2.5});
    check_voigt<1>(output, {2.5});
  }

  SECTION("2D plane stress Tensor")
  {
    dealii::Tensor<2, 2, double> strain;
    strain[0][0] = 1.0;
    strain[1][1] = 2.0;
    strain[0][1] = 3.0;
    strain[1][0] = 5.0;

    const auto returned = Mechanics::strain_to_voigt<2, plane_stress>(strain);

    dealii::Tensor<1, 3, double> output;
    Mechanics::strain_to_voigt<2, plane_stress>(strain, output);

    // Engineering shear strain: gamma_xy = epsilon_xy + epsilon_yx.
    check_voigt<3>(returned, {1.0, 2.0, 8.0});
    check_voigt<3>(output, {1.0, 2.0, 8.0});
  }

  SECTION("2D plane stress SymmetricTensor")
  {
    dealii::SymmetricTensor<2, 2, double> strain;
    strain[0][0] = 1.0;
    strain[1][1] = 2.0;
    strain[0][1] = 4.0;

    const auto returned = Mechanics::strain_to_voigt<2, plane_stress>(strain);

    dealii::Tensor<1, 3, double> output;
    Mechanics::strain_to_voigt<2, plane_stress>(strain, output);

    check_voigt<3>(returned, {1.0, 2.0, 8.0});
    check_voigt<3>(output, {1.0, 2.0, 8.0});
  }

  SECTION("2D plane strain Tensor")
  {
    dealii::Tensor<2, 2, double> strain;
    strain[0][0] = 1.0;
    strain[1][1] = 2.0;
    strain[0][1] = 3.0;
    strain[1][0] = 5.0;

    const double strain_zz = 6.0;

    const auto returned = Mechanics::strain_to_voigt<2, plane_strain>(strain, strain_zz);

    dealii::Tensor<1, 4, double> output;
    Mechanics::strain_to_voigt<2, plane_strain>(strain, strain_zz, output);

    check_voigt<4>(returned, {1.0, 2.0, 6.0, 8.0});
    check_voigt<4>(output, {1.0, 2.0, 6.0, 8.0});
  }

  SECTION("2D plane strain SymmetricTensor")
  {
    dealii::SymmetricTensor<2, 2, double> strain;
    strain[0][0] = 1.0;
    strain[1][1] = 2.0;
    strain[0][1] = 4.0;

    const double strain_zz = 6.0;

    const auto returned = Mechanics::strain_to_voigt<2, plane_strain>(strain, strain_zz);

    dealii::Tensor<1, 4, double> output;
    Mechanics::strain_to_voigt<2, plane_strain>(strain, strain_zz, output);

    check_voigt<4>(returned, {1.0, 2.0, 6.0, 8.0});
    check_voigt<4>(output, {1.0, 2.0, 6.0, 8.0});
  }

  SECTION("3D Tensor")
  {
    dealii::Tensor<2, 3, double> strain;

    strain[0][0] = 1.0;
    strain[1][1] = 2.0;
    strain[2][2] = 3.0;

    strain[1][2] = 4.0;
    strain[2][1] = 5.0;
    strain[0][2] = 6.0;
    strain[2][0] = 7.0;
    strain[0][1] = 8.0;
    strain[1][0] = 9.0;

    const auto returned = Mechanics::strain_to_voigt<3>(strain);

    dealii::Tensor<1, 6, double> output;
    Mechanics::strain_to_voigt<3>(strain, output);

    check_voigt<6>(returned, {1.0, 2.0, 3.0, 9.0, 13.0, 17.0});
    check_voigt<6>(output, {1.0, 2.0, 3.0, 9.0, 13.0, 17.0});
  }

  SECTION("3D SymmetricTensor")
  {
    dealii::SymmetricTensor<2, 3, double> strain;

    strain[0][0] = 1.0;
    strain[1][1] = 2.0;
    strain[2][2] = 3.0;
    strain[1][2] = 4.0;
    strain[0][2] = 5.0;
    strain[0][1] = 6.0;

    const auto returned = Mechanics::strain_to_voigt<3>(strain);

    dealii::Tensor<1, 6, double> output;
    Mechanics::strain_to_voigt<3>(strain, output);

    check_voigt<6>(returned, {1.0, 2.0, 3.0, 8.0, 10.0, 12.0});
    check_voigt<6>(output, {1.0, 2.0, 3.0, 8.0, 10.0, 12.0});
  }
}

TEST_CASE("Mechanics voigt_to_strain", "[mechanics][strain][voigt]")
{
  SECTION("1D")
  {
    const auto voigt = make_voigt<1>({2.5});

    dealii::Tensor<2, 1, double> tensor;
    Mechanics::voigt_to_strain<1>(voigt, tensor);

    dealii::SymmetricTensor<2, 1, double> symmetric_tensor;
    Mechanics::voigt_to_strain<1>(voigt, symmetric_tensor);

    const auto returned = Mechanics::voigt_to_strain<1>(voigt);

    check_close(tensor[0][0], 2.5);
    check_close(symmetric_tensor[0][0], 2.5);
    check_close(returned[0][0], 2.5);
  }

  SECTION("2D plane stress")
  {
    const auto voigt = make_voigt<3>({1.0, 2.0, 8.0});

    dealii::Tensor<2, 2, double> tensor;
    Mechanics::voigt_to_strain<2, plane_stress>(voigt, tensor);

    dealii::SymmetricTensor<2, 2, double> symmetric_tensor;
    Mechanics::voigt_to_strain<2, plane_stress>(voigt, symmetric_tensor);

    const auto returned = Mechanics::voigt_to_strain<2, plane_stress>(voigt);

    check_close(tensor[0][0], 1.0);
    check_close(tensor[1][1], 2.0);
    check_close(tensor[0][1], 4.0);
    check_close(tensor[1][0], 4.0);

    check_rank_2_tensor<2>(symmetric_tensor, tensor);
    check_rank_2_tensor<2>(returned, tensor);
  }

  SECTION("2D plane strain")
  {
    const auto voigt = make_voigt<4>({1.0, 2.0, 6.0, 8.0});

    dealii::Tensor<2, 2, double> tensor;
    double                       component_zz = 0.0;

    Mechanics::voigt_to_strain<2, plane_strain>(voigt, tensor, component_zz);

    dealii::SymmetricTensor<2, 2, double> symmetric_tensor;
    double                                symmetric_component_zz = 0.0;

    Mechanics::voigt_to_strain<2, plane_strain>(voigt,
                                                symmetric_tensor,
                                                symmetric_component_zz);

    const auto [returned_tensor, returned_component_zz] =
      Mechanics::voigt_to_strain<2, plane_strain>(voigt);

    check_close(tensor[0][0], 1.0);
    check_close(tensor[1][1], 2.0);
    check_close(tensor[0][1], 4.0);
    check_close(tensor[1][0], 4.0);
    check_close(component_zz, 6.0);

    check_rank_2_tensor<2>(symmetric_tensor, tensor);
    check_rank_2_tensor<2>(returned_tensor, tensor);
    check_close(symmetric_component_zz, 6.0);
    check_close(returned_component_zz, 6.0);
  }

  SECTION("3D")
  {
    const auto voigt = make_voigt<6>({1.0, 2.0, 3.0, 8.0, 10.0, 12.0});

    dealii::Tensor<2, 3, double> tensor;
    Mechanics::voigt_to_strain<3>(voigt, tensor);

    dealii::SymmetricTensor<2, 3, double> symmetric_tensor;
    Mechanics::voigt_to_strain<3>(voigt, symmetric_tensor);

    const auto returned = Mechanics::voigt_to_strain<3>(voigt);

    check_close(tensor[0][0], 1.0);
    check_close(tensor[1][1], 2.0);
    check_close(tensor[2][2], 3.0);

    check_close(tensor[1][2], 4.0);
    check_close(tensor[2][1], 4.0);
    check_close(tensor[0][2], 5.0);
    check_close(tensor[2][0], 5.0);
    check_close(tensor[0][1], 6.0);
    check_close(tensor[1][0], 6.0);

    check_rank_2_tensor<3>(symmetric_tensor, tensor);
    check_rank_2_tensor<3>(returned, tensor);
  }
}

TEST_CASE("Mechanics strain Voigt round trips", "[mechanics][strain][voigt][round-trip]")
{
  SECTION("2D plane stress")
  {
    dealii::SymmetricTensor<2, 2, double> original;
    original[0][0] = 1.0;
    original[1][1] = 2.0;
    original[0][1] = 3.0;

    const auto voigt = Mechanics::strain_to_voigt<2, plane_stress>(original);

    const auto recovered = Mechanics::voigt_to_strain<2, plane_stress>(voigt);

    check_rank_2_tensor<2>(recovered, original);
  }

  SECTION("2D plane strain")
  {
    dealii::SymmetricTensor<2, 2, double> original;
    original[0][0] = 1.0;
    original[1][1] = 2.0;
    original[0][1] = 3.0;

    const double original_zz = 4.0;

    const auto voigt = Mechanics::strain_to_voigt<2, plane_strain>(original, original_zz);

    const auto [recovered, recovered_zz] =
      Mechanics::voigt_to_strain<2, plane_strain>(voigt);

    check_rank_2_tensor<2>(recovered, original);
    check_close(recovered_zz, original_zz);
  }

  SECTION("3D")
  {
    dealii::SymmetricTensor<2, 3, double> original;
    original[0][0] = 1.0;
    original[1][1] = 2.0;
    original[2][2] = 3.0;
    original[1][2] = 4.0;
    original[0][2] = 5.0;
    original[0][1] = 6.0;

    const auto voigt     = Mechanics::strain_to_voigt<3>(original);
    const auto recovered = Mechanics::voigt_to_strain<3>(voigt);

    check_rank_2_tensor<3>(recovered, original);
  }
}

TEST_CASE("Mechanics stress_to_voigt", "[mechanics][stress][voigt]")
{
  SECTION("1D Tensor")
  {
    dealii::Tensor<2, 1, double> stress;
    stress[0][0] = 2.5;

    const auto returned = Mechanics::stress_to_voigt<1>(stress);

    dealii::Tensor<1, 1, double> output;
    Mechanics::stress_to_voigt<1>(stress, output);

    check_voigt<1>(returned, {2.5});
    check_voigt<1>(output, {2.5});
  }

  SECTION("1D SymmetricTensor")
  {
    dealii::SymmetricTensor<2, 1, double> stress;
    stress[0][0] = 2.5;

    const auto returned = Mechanics::stress_to_voigt<1>(stress);

    dealii::Tensor<1, 1, double> output;
    Mechanics::stress_to_voigt<1>(stress, output);

    check_voigt<1>(returned, {2.5});
    check_voigt<1>(output, {2.5});
  }

  SECTION("2D plane stress Tensor")
  {
    dealii::Tensor<2, 2, double> stress;
    stress[0][0] = 1.0;
    stress[1][1] = 2.0;
    stress[0][1] = 3.0;
    stress[1][0] = 5.0;

    const auto returned = Mechanics::stress_to_voigt<2, plane_stress>(stress);

    dealii::Tensor<1, 3, double> output;
    Mechanics::stress_to_voigt<2, plane_stress>(stress, output);

    check_voigt<3>(returned, {1.0, 2.0, 4.0});
    check_voigt<3>(output, {1.0, 2.0, 4.0});
  }

  SECTION("2D plane stress SymmetricTensor")
  {
    dealii::SymmetricTensor<2, 2, double> stress;
    stress[0][0] = 1.0;
    stress[1][1] = 2.0;
    stress[0][1] = 4.0;

    const auto returned = Mechanics::stress_to_voigt<2, plane_stress>(stress);

    dealii::Tensor<1, 3, double> output;
    Mechanics::stress_to_voigt<2, plane_stress>(stress, output);

    check_voigt<3>(returned, {1.0, 2.0, 4.0});
    check_voigt<3>(output, {1.0, 2.0, 4.0});
  }

  SECTION("2D plane strain Tensor")
  {
    dealii::Tensor<2, 2, double> stress;
    stress[0][0] = 1.0;
    stress[1][1] = 2.0;
    stress[0][1] = 4.0;
    stress[1][0] = 4.0;

    const double stress_zz = 6.0;

    const auto returned = Mechanics::stress_to_voigt<2, plane_strain>(stress, stress_zz);

    dealii::Tensor<1, 4, double> output;
    Mechanics::stress_to_voigt<2, plane_strain>(stress, stress_zz, output);

    check_voigt<4>(returned, {1.0, 2.0, 6.0, 4.0});
    check_voigt<4>(output, {1.0, 2.0, 6.0, 4.0});
  }

  SECTION("2D plane strain SymmetricTensor")
  {
    dealii::SymmetricTensor<2, 2, double> stress;
    stress[0][0] = 1.0;
    stress[1][1] = 2.0;
    stress[0][1] = 4.0;

    const double stress_zz = 6.0;

    const auto returned = Mechanics::stress_to_voigt<2, plane_strain>(stress, stress_zz);

    dealii::Tensor<1, 4, double> output;
    Mechanics::stress_to_voigt<2, plane_strain>(stress, stress_zz, output);

    check_voigt<4>(returned, {1.0, 2.0, 6.0, 4.0});
    check_voigt<4>(output, {1.0, 2.0, 6.0, 4.0});
  }

  SECTION("3D Tensor")
  {
    dealii::Tensor<2, 3, double> stress;

    stress[0][0] = 1.0;
    stress[1][1] = 2.0;
    stress[2][2] = 3.0;

    stress[1][2] = stress[2][1] = 4.0;
    stress[0][2] = stress[2][0] = 5.0;
    stress[0][1] = stress[1][0] = 6.0;

    const auto returned = Mechanics::stress_to_voigt<3>(stress);

    dealii::Tensor<1, 6, double> output;
    Mechanics::stress_to_voigt<3>(stress, output);

    check_voigt<6>(returned, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    check_voigt<6>(output, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
  }

  SECTION("3D SymmetricTensor")
  {
    dealii::SymmetricTensor<2, 3, double> stress;

    stress[0][0] = 1.0;
    stress[1][1] = 2.0;
    stress[2][2] = 3.0;
    stress[1][2] = 4.0;
    stress[0][2] = 5.0;
    stress[0][1] = 6.0;

    const auto returned = Mechanics::stress_to_voigt<3>(stress);

    dealii::Tensor<1, 6, double> output;
    Mechanics::stress_to_voigt<3>(stress, output);

    check_voigt<6>(returned, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    check_voigt<6>(output, {1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
  }
}

TEST_CASE("Mechanics voigt_to_stress", "[mechanics][stress][voigt]")
{
  SECTION("1D")
  {
    const auto voigt = make_voigt<1>({2.5});

    dealii::Tensor<2, 1, double> tensor;
    Mechanics::voigt_to_stress<1>(voigt, tensor);

    dealii::SymmetricTensor<2, 1, double> symmetric_tensor;
    Mechanics::voigt_to_stress<1>(voigt, symmetric_tensor);

    const auto returned = Mechanics::voigt_to_stress<1>(voigt);

    check_close(tensor[0][0], 2.5);
    check_close(symmetric_tensor[0][0], 2.5);
    check_close(returned[0][0], 2.5);
  }

  SECTION("2D plane stress")
  {
    const auto voigt = make_voigt<3>({1.0, 2.0, 4.0});

    dealii::Tensor<2, 2, double> tensor;
    Mechanics::voigt_to_stress<2, plane_stress>(voigt, tensor);

    dealii::SymmetricTensor<2, 2, double> symmetric_tensor;
    Mechanics::voigt_to_stress<2, plane_stress>(voigt, symmetric_tensor);

    const auto returned = Mechanics::voigt_to_stress<2, plane_stress>(voigt);

    check_close(tensor[0][0], 1.0);
    check_close(tensor[1][1], 2.0);
    check_close(tensor[0][1], 4.0);
    check_close(tensor[1][0], 4.0);

    check_rank_2_tensor<2>(symmetric_tensor, tensor);
    check_rank_2_tensor<2>(returned, tensor);
  }

  SECTION("2D plane strain")
  {
    const auto voigt = make_voigt<4>({1.0, 2.0, 6.0, 4.0});

    dealii::Tensor<2, 2, double> tensor;
    double                       component_zz = 0.0;

    Mechanics::voigt_to_stress<2, plane_strain>(voigt, tensor, component_zz);

    dealii::SymmetricTensor<2, 2, double> symmetric_tensor;
    double                                symmetric_component_zz = 0.0;

    Mechanics::voigt_to_stress<2, plane_strain>(voigt,
                                                symmetric_tensor,
                                                symmetric_component_zz);

    const auto [returned_tensor, returned_component_zz] =
      Mechanics::voigt_to_stress<2, plane_strain>(voigt);

    check_close(tensor[0][0], 1.0);
    check_close(tensor[1][1], 2.0);
    check_close(tensor[0][1], 4.0);
    check_close(tensor[1][0], 4.0);
    check_close(component_zz, 6.0);

    check_rank_2_tensor<2>(symmetric_tensor, tensor);
    check_rank_2_tensor<2>(returned_tensor, tensor);
    check_close(symmetric_component_zz, 6.0);
    check_close(returned_component_zz, 6.0);
  }

  SECTION("3D")
  {
    const auto voigt = make_voigt<6>({1.0, 2.0, 3.0, 4.0, 5.0, 6.0});

    dealii::Tensor<2, 3, double> tensor;
    Mechanics::voigt_to_stress<3>(voigt, tensor);

    dealii::SymmetricTensor<2, 3, double> symmetric_tensor;
    Mechanics::voigt_to_stress<3>(voigt, symmetric_tensor);

    const auto returned = Mechanics::voigt_to_stress<3>(voigt);

    check_close(tensor[0][0], 1.0);
    check_close(tensor[1][1], 2.0);
    check_close(tensor[2][2], 3.0);

    check_close(tensor[1][2], 4.0);
    check_close(tensor[2][1], 4.0);
    check_close(tensor[0][2], 5.0);
    check_close(tensor[2][0], 5.0);
    check_close(tensor[0][1], 6.0);
    check_close(tensor[1][0], 6.0);

    check_rank_2_tensor<3>(symmetric_tensor, tensor);
    check_rank_2_tensor<3>(returned, tensor);
  }
}

TEST_CASE("Mechanics stress Voigt round trips", "[mechanics][stress][voigt][round-trip]")
{
  SECTION("Stress tensor to Voigt to stress tensor")
  {
    dealii::SymmetricTensor<2, 2, double> original;
    original[0][0] = 1.0;
    original[1][1] = 2.0;
    original[0][1] = 3.0;

    const auto voigt = Mechanics::stress_to_voigt<2, plane_stress>(original);

    const auto recovered = Mechanics::voigt_to_stress<2, plane_stress>(voigt);

    check_rank_2_tensor<2>(recovered, original);
  }

  SECTION("Voigt to stress tensor to Voigt")
  {
    const auto original = make_voigt<3>({1.0, 2.0, 3.0});

    const auto tensor = Mechanics::voigt_to_stress<2, plane_stress>(original);

    const auto recovered = Mechanics::stress_to_voigt<2, plane_stress>(tensor);

    check_voigt<3>(recovered, {1.0, 2.0, 3.0});
  }

  SECTION("Plane strain")
  {
    dealii::SymmetricTensor<2, 2, double> original;
    original[0][0] = 1.0;
    original[1][1] = 2.0;
    original[0][1] = 3.0;

    const double original_zz = 4.0;

    const auto voigt = Mechanics::stress_to_voigt<2, plane_strain>(original, original_zz);

    const auto [recovered, recovered_zz] =
      Mechanics::voigt_to_stress<2, plane_strain>(voigt);

    check_rank_2_tensor<2>(recovered, original);
    check_close(recovered_zz, original_zz);
  }
}

TEST_CASE("Mechanics isotropic stiffness", "[mechanics][stiffness][isotropic]")
{
  SECTION("1D")
  {
    const auto stiffness = Mechanics::stiffness_isotropic<1>(1.0, 0.25);

    check_close(stiffness[0][0], 1.0);
  }

  SECTION("2D plane stress")
  {
    const auto stiffness = Mechanics::stiffness_isotropic<2, plane_stress>(1.0, 0.25);

    dealii::Tensor<2, 3, double> expected;

    expected[0][0] = 16.0 / 15.0;
    expected[1][1] = 16.0 / 15.0;
    expected[0][1] = 4.0 / 15.0;
    expected[1][0] = 4.0 / 15.0;
    expected[2][2] = 2.0 / 5.0;

    check_matrix<3>(stiffness, expected);
  }

  SECTION("2D plane strain")
  {
    const auto stiffness = Mechanics::stiffness_isotropic<2, plane_strain>(1.0, 0.25);

    dealii::Tensor<2, 4, double> expected;

    expected[0][0] = 1.2;
    expected[1][1] = 1.2;
    expected[2][2] = 1.2;

    expected[0][1] = expected[1][0] = 0.4;
    expected[0][2] = expected[2][0] = 0.4;
    expected[1][2] = expected[2][1] = 0.4;

    expected[3][3] = 0.4;

    check_matrix<4>(stiffness, expected);
  }

  SECTION("3D")
  {
    const auto stiffness = Mechanics::stiffness_isotropic<3>(1.0, 0.25);

    dealii::Tensor<2, 6, double> expected;

    expected[0][0] = 1.2;
    expected[1][1] = 1.2;
    expected[2][2] = 1.2;

    expected[0][1] = expected[1][0] = 0.4;
    expected[0][2] = expected[2][0] = 0.4;
    expected[1][2] = expected[2][1] = 0.4;

    expected[3][3] = 0.4;
    expected[4][4] = 0.4;
    expected[5][5] = 0.4;

    check_matrix<6>(stiffness, expected);
  }
}

#ifdef DEBUG
TEST_CASE("Mechanics isotropic stiffness rejects invalid constants",
          "[mechanics][stiffness][isotropic][exceptions]")
{
  CHECK_THROWS(Mechanics::stiffness_isotropic<1>(0.0, 0.25));
  CHECK_THROWS(Mechanics::stiffness_isotropic<1>(-1.0, 0.25));

  CHECK_THROWS(Mechanics::stiffness_isotropic<1>(1.0, -1.0));
  CHECK_THROWS(Mechanics::stiffness_isotropic<1>(1.0, 0.5));
}
#endif

TEST_CASE("Mechanics orthotropic stiffness", "[mechanics][stiffness][orthotropic]")
{
  SECTION("2D plane stress isotropic limit")
  {
    const auto orthotropic =
      Mechanics::stiffness_orthotropic<2, plane_stress>(1.0, 1.0, 0.25, 0.4);

    const auto isotropic = Mechanics::stiffness_isotropic<2, plane_stress>(1.0, 0.25);

    check_matrix<3>(orthotropic, isotropic);
  }

  SECTION("2D plane strain isotropic limit")
  {
    const auto orthotropic = Mechanics::stiffness_orthotropic<2, plane_strain>(1.0,
                                                                               1.0,
                                                                               1.0,
                                                                               0.25,
                                                                               0.25,
                                                                               0.25,
                                                                               0.4);

    const auto isotropic = Mechanics::stiffness_isotropic<2, plane_strain>(1.0, 0.25);

    check_matrix<4>(orthotropic, isotropic);
  }

  SECTION("3D isotropic limit")
  {
    const auto orthotropic = Mechanics::stiffness_orthotropic<3, three_dimensional>(1.0,
                                                                                    1.0,
                                                                                    1.0,
                                                                                    0.25,
                                                                                    0.25,
                                                                                    0.25,
                                                                                    0.4,
                                                                                    0.4,
                                                                                    0.4);

    const auto isotropic = Mechanics::stiffness_isotropic<3>(1.0, 0.25);

    check_matrix<6>(orthotropic, isotropic);
  }
}

#ifdef DEBUG
TEST_CASE("Mechanics orthotropic stiffness rejects invalid constants",
          "[mechanics][stiffness][orthotropic][exceptions]")
{
  CHECK_THROWS((Mechanics::stiffness_orthotropic<2, plane_stress>(0.0, 1.0, 0.25, 0.4)));

  CHECK_THROWS((Mechanics::stiffness_orthotropic<2, plane_stress>(1.0, 0.0, 0.25, 0.4)));

  CHECK_THROWS((Mechanics::stiffness_orthotropic<2, plane_stress>(1.0, 1.0, 0.25, 0.0)));

  CHECK_THROWS((Mechanics::stiffness_orthotropic<2, plane_stress>(1.0, 1.0, 1.0, 0.4)));
}
#endif

TEST_CASE("Mechanics extract_plane_strain_stiffness",
          "[mechanics][stiffness][plane-strain]")
{
  const auto stiffness_3d = Mechanics::stiffness_isotropic<3>(1.0, 0.25);

  const auto extracted = Mechanics::extract_plane_strain_stiffness(stiffness_3d);

  const auto expected = Mechanics::stiffness_isotropic<2, plane_strain>(1.0, 0.25);

  check_matrix<4>(extracted, expected);
}

TEST_CASE("Mechanics compute_stress for plane stress",
          "[mechanics][compute-stress][plane-stress]")
{
  const auto stiffness = Mechanics::stiffness_isotropic<2, plane_stress>(1.0, 0.25);

  SECTION("Voigt strain and stress")
  {
    const auto strain = make_voigt<3>({1.0, 0.0, 0.0});

    dealii::Tensor<1, 3, double> stress;

    Mechanics::compute_stress<2, plane_stress>(stiffness, strain, stress);

    check_voigt<3>(stress, {16.0 / 15.0, 4.0 / 15.0, 0.0});
  }

  SECTION("Tensor strain and stress")
  {
    dealii::Tensor<2, 2, double> strain;
    strain[0][0] = 1.0;

    dealii::Tensor<2, 2, double> stress;

    Mechanics::compute_stress<2, plane_stress>(stiffness, strain, stress);

    check_close(stress[0][0], 16.0 / 15.0);
    check_close(stress[1][1], 4.0 / 15.0);
    check_close(stress[0][1], 0.0);
    check_close(stress[1][0], 0.0);
  }

  SECTION("Tensor shear strain")
  {
    dealii::Tensor<2, 2, double> strain;
    strain[0][1] = 0.5;
    strain[1][0] = 0.5;

    dealii::Tensor<2, 2, double> stress;

    Mechanics::compute_stress<2, plane_stress>(stiffness, strain, stress);

    check_close(stress[0][0], 0.0);
    check_close(stress[1][1], 0.0);
    check_close(stress[0][1], 0.4);
    check_close(stress[1][0], 0.4);
  }
}

TEST_CASE("Mechanics compute_stress for plane strain",
          "[mechanics][compute-stress][plane-strain]")
{
  const auto stiffness = Mechanics::stiffness_isotropic<2, plane_strain>(1.0, 0.25);

  SECTION("Voigt strain and stress")
  {
    const auto strain = make_voigt<4>({1.0, 0.0, 0.0, 0.0});

    dealii::Tensor<1, 4, double> stress;

    Mechanics::compute_stress<2, plane_strain>(stiffness, strain, stress);

    check_voigt<4>(stress, {1.2, 0.4, 0.4, 0.0});
  }

  SECTION("Tensor strain and tensor stress")
  {
    dealii::Tensor<2, 2, double> strain;
    strain[0][0] = 1.0;

    constexpr double strain_zz = 0.0;

    dealii::Tensor<2, 2, double> stress;
    double                       stress_zz = 0.0;

    Mechanics::compute_stress<2, plane_strain>(stiffness,
                                               strain,
                                               strain_zz,
                                               stress,
                                               stress_zz);

    check_close(stress[0][0], 1.2);
    check_close(stress[1][1], 0.4);
    check_close(stress[0][1], 0.0);
    check_close(stress[1][0], 0.0);
    check_close(stress_zz, 0.4);
  }

  SECTION("Tensor strain and Voigt stress")
  {
    dealii::Tensor<2, 2, double> strain;
    strain[0][0] = 1.0;

    constexpr double strain_zz = 0.0;

    dealii::Tensor<1, 4, double> stress;
    dealii::Tensor<1, 4, double> stress2;

    Mechanics::compute_stress<2, plane_strain>(stiffness, strain, strain_zz, stress);
    Mechanics::compute_stress<2, plane_strain>(stiffness, strain, stress2);

    check_voigt<4>(stress, {1.2, 0.4, 0.4, 0.0});
    check_voigt<4>(stress2, {1.2, 0.4, 0.4, 0.0});
  }

  SECTION("Nonzero out-of-plane strain")
  {
    dealii::Tensor<2, 2, double> strain;

    constexpr double strain_zz = 1.0;

    dealii::Tensor<2, 2, double> stress;
    double                       stress_zz = 0.0;

    Mechanics::compute_stress<2, plane_strain>(stiffness,
                                               strain,
                                               strain_zz,
                                               stress,
                                               stress_zz);

    check_close(stress[0][0], 0.4);
    check_close(stress[1][1], 0.4);
    check_close(stress[0][1], 0.0);
    check_close(stress[1][0], 0.0);
    check_close(stress_zz, 1.2);
  }
}

TEST_CASE("Mechanics compute_stress in 1D and 3D", "[mechanics][compute-stress]")
{
  SECTION("1D")
  {
    const auto stiffness = Mechanics::stiffness_isotropic<1>(2.0, 0.25);

    const auto strain = make_voigt<1>({3.0});

    dealii::Tensor<1, 1, double> stress;

    Mechanics::compute_stress<1, three_dimensional>(stiffness, strain, stress);

    check_voigt<1>(stress, {6.0});
  }

  SECTION("3D uniaxial strain")
  {
    const auto stiffness = Mechanics::stiffness_isotropic<3>(1.0, 0.25);

    const auto strain = make_voigt<6>({1.0, 0.0, 0.0, 0.0, 0.0, 0.0});

    dealii::Tensor<1, 6, double> stress;

    Mechanics::compute_stress<3, three_dimensional>(stiffness, strain, stress);

    check_voigt<6>(stress, {1.2, 0.4, 0.4, 0.0, 0.0, 0.0});
  }

  SECTION("3D tensor shear strain")
  {
    const auto stiffness = Mechanics::stiffness_isotropic<3>(1.0, 0.25);

    dealii::Tensor<2, 3, double> strain;
    strain[0][1] = 0.5;
    strain[1][0] = 0.5;

    dealii::Tensor<2, 3, double> stress;

    Mechanics::compute_stress<3, three_dimensional>(stiffness, strain, stress);

    check_close(stress[0][0], 0.0);
    check_close(stress[1][1], 0.0);
    check_close(stress[2][2], 0.0);
    check_close(stress[0][1], 0.4);
    check_close(stress[1][0], 0.4);
  }
}

TEST_CASE("Mechanics strain energy", "[mechanics][strain-energy]")
{
  SECTION("1D")
  {
    const auto stress   = make_voigt<1>({4.0});
    const auto strain_e = make_voigt<1>({3.0});

    const double energy =
      Mechanics::strain_energy<1, three_dimensional>(stress, strain_e);

    check_close(energy, 6.0);
  }

  SECTION("2D plane stress")
  {
    const auto stress   = make_voigt<3>({2.0, 4.0, 6.0});
    const auto strain_e = make_voigt<3>({1.0, 0.5, 0.25});

    const double energy = Mechanics::strain_energy<2, plane_stress>(stress, strain_e);

    check_close(energy, 2.75);
  }

  SECTION("2D plane strain")
  {
    const auto stress   = make_voigt<4>({2.0, 4.0, 6.0, 8.0});
    const auto strain_e = make_voigt<4>({1.0, 0.5, 0.25, 0.125});

    const double energy = Mechanics::strain_energy<2, plane_strain>(stress, strain_e);

    check_close(energy, 3.25);
  }

  SECTION("3D")
  {
    const auto stress   = make_voigt<6>({1.0, 2.0, 3.0, 4.0, 5.0, 6.0});
    const auto strain_e = make_voigt<6>({6.0, 5.0, 4.0, 3.0, 2.0, 1.0});

    const double energy =
      Mechanics::strain_energy<3, three_dimensional>(stress, strain_e);

    check_close(energy, 28.0);
  }
}

TEST_CASE("Mechanics von Mises stress", "[mechanics][stress-mises]")
{
  SECTION("1D returns the magnitude of the stress")
  {
    const auto stress = make_voigt<1>({-7.0});

    const double mises = Mechanics::stress_mises<1, three_dimensional>(stress);

    check_close(mises, 7.0);
  }

  SECTION("2D plane stress")
  {
    const auto stress = make_voigt<3>({100.0, 40.0, 30.0});

    const double mises = Mechanics::stress_mises<2, plane_stress>(stress);

    check_close(mises, std::sqrt(10300.0));
  }

  SECTION("2D plane strain includes sigma zz")
  {
    const auto stress = make_voigt<4>({100.0, 40.0, 20.0, 30.0});

    const double mises = Mechanics::stress_mises<2, plane_strain>(stress);

    check_close(mises, std::sqrt(7900.0));
  }

  SECTION("3D")
  {
    const auto stress = make_voigt<6>({100.0, 40.0, 20.0, 10.0, 20.0, 30.0});

    const double mises = Mechanics::stress_mises<3, three_dimensional>(stress);

    check_close(mises, std::sqrt(9400.0));
  }

  SECTION("3D hydrostatic stress has zero von Mises stress")
  {
    const auto stress = make_voigt<6>({25.0, 25.0, 25.0, 0.0, 0.0, 0.0});

    const double mises = Mechanics::stress_mises<3, three_dimensional>(stress);

    check_close(mises, 0.0);
  }

  SECTION("3D pure shear")
  {
    const auto stress = make_voigt<6>({0.0, 0.0, 0.0, 0.0, 0.0, 8.0});

    const double mises = Mechanics::stress_mises<3, three_dimensional>(stress);

    check_close(mises, 8.0 * std::numbers::sqrt3);
  }
}

TEST_CASE("Mechanics principal stress", "[mechanics][stress-principal]")
{
  SECTION("1D")
  {
    const auto stress = make_voigt<1>({-12.0});

    const auto principal = Mechanics::stress_principal<1, three_dimensional>(stress);

    check_voigt<1>(principal, {-12.0});
  }

  SECTION("2D plane stress")
  {
    // The in-plane stress tensor is
    // [ 5  2 ]
    // [ 2  2 ]
    // Its eigenvalues are 6 and 1.
    const auto stress = make_voigt<3>({5.0, 2.0, 2.0});

    const auto principal = Mechanics::stress_principal<2, plane_stress>(stress);

    check_voigt<2>(principal, {6.0, 1.0});
  }

  SECTION("2D plane strain uses the in-plane stress tensor")
  {
    // sigma_zz is intentionally large. The function should return the
    // principal stresses of the in-plane 2x2 tensor and ignore sigma_zz.
    const auto stress = make_voigt<4>({5.0, 2.0, 123.0, 2.0});

    const auto principal = Mechanics::stress_principal<2, plane_strain>(stress);

    check_voigt<2>(principal, {6.0, 1.0});
  }

  SECTION("2D hydrostatic in-plane stress")
  {
    const auto stress = make_voigt<3>({9.0, 9.0, 0.0});

    const auto principal = Mechanics::stress_principal<2, plane_stress>(stress);

    check_voigt<2>(principal, {9.0, 9.0});
  }
}
