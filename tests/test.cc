#include "../src/models/mechanics/computeStress.h"

#define CATCH_CONFIG_MAIN
#include "../contrib/catch/catch.hpp"

TEST_CASE("Compute Stress")
{
  const double tolerance = 1.0e-6;

  SECTION("1D")
  {
    const unsigned int dim = 1;

    dealii::Tensor<2, 2 * dim - 1 + dim / 3, dealii::VectorizedArray<double>> CIJ;
    dealii::VectorizedArray<double> R[dim][dim], ux[dim][dim];

    CIJ[0][0] = 2.5;

    ux[0][0] = 1.0;

    computeStress<dim>(CIJ, ux, R);

    REQUIRE(R[0][0][0] == Approx(2.5).epsilon(tolerance));
  }
  SECTION("2D")
  {
    const unsigned int dim = 2;

    dealii::Tensor<2, 2 * dim - 1 + dim / 3, dealii::VectorizedArray<double>> CIJ;
    dealii::VectorizedArray<double> R[dim][dim], ux[dim][dim];

    CIJ[0][0] = 6.8;
    CIJ[1][0] = 2.5;
    CIJ[2][0] = 4.0;
    CIJ[0][1] = 1.2;
    CIJ[1][1] = 10.1;
    CIJ[2][1] = 3.7;
    CIJ[0][2] = 9.3;
    CIJ[1][2] = 2.7;
    CIJ[2][2] = 8.8;

    ux[0][0] = 1.0;
    ux[1][0] = 2.0;
    ux[0][1] = 3.0;
    ux[1][1] = 4.0;

    computeStress<dim>(CIJ, ux, R);

    REQUIRE(R[0][0][0] == Approx(58.1).epsilon(tolerance));
    REQUIRE(R[1][1][0] == Approx(56.4).epsilon(tolerance));
    REQUIRE(R[0][1][0] == Approx(62.8).epsilon(tolerance));
    REQUIRE(R[1][0][0] == Approx(62.8).epsilon(tolerance));
  }
  SECTION("3D")
  {
    const unsigned int dim = 3;

    dealii::Tensor<2, 2 * dim - 1 + dim / 3, dealii::VectorizedArray<double>> CIJ;
    dealii::VectorizedArray<double> R[dim][dim], ux[dim][dim];

    CIJ[0][0] = 1.1;
    CIJ[1][1] = 7.7;
    CIJ[2][2] = 6.6;
    CIJ[3][3] = 3.3;
    CIJ[4][4] = 11.6;
    CIJ[5][5] = 19.5;
    CIJ[0][1] = CIJ[1][0] = 9.5;
    CIJ[0][2] = CIJ[2][0] = 2.1;
    CIJ[0][3] = CIJ[3][0] = 1.5;
    CIJ[0][4] = CIJ[4][0] = 9.2;
    CIJ[0][5] = CIJ[5][0] = 18.6;
    CIJ[1][2] = CIJ[2][1] = 5.6;
    CIJ[1][3] = CIJ[3][1] = 4.7;
    CIJ[1][4] = CIJ[4][1] = 6.4;
    CIJ[1][5] = CIJ[5][1] = 5.9;
    CIJ[2][3] = CIJ[3][2] = 15.5;
    CIJ[2][4] = CIJ[4][2] = 63.1;
    CIJ[2][5] = CIJ[5][2] = 50.0;
    CIJ[3][4] = CIJ[4][3] = 92.5;
    CIJ[3][5] = CIJ[5][3] = 1.3;
    CIJ[4][5] = CIJ[5][4] = 23.2;

    ux[0][0] = 1.0;
    ux[1][0] = 2.0;
    ux[2][0] = 3.0;
    ux[0][1] = 4.0;
    ux[1][1] = 5.0;
    ux[2][1] = 6.0;
    ux[0][2] = 7.0;
    ux[1][2] = 8.0;
    ux[2][2] = 9.0;

    computeStress<dim>(CIJ, ux, R);

    REQUIRE(R[0][0][0] == Approx(292.1).epsilon(tolerance));
    REQUIRE(R[1][1][0] == Approx(263.6).epsilon(tolerance));
    REQUIRE(R[2][2][0] == Approx(1237.5).epsilon(tolerance));
    REQUIRE(R[1][2][0] == Approx(1143.5).epsilon(tolerance));
    REQUIRE(R[2][1][0] == Approx(1143.5).epsilon(tolerance));
    REQUIRE(R[0][2][0] == Approx(2159.3).epsilon(tolerance));
    REQUIRE(R[2][0][0] == Approx(2159.3).epsilon(tolerance));
    REQUIRE(R[0][1][0] == Approx(865.3).epsilon(tolerance));
    REQUIRE(R[1][0][0] == Approx(865.3).epsilon(tolerance));
  }
}