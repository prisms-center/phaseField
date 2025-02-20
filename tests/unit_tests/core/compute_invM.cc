// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/grid/grid_generator.h>

#include "catch.hpp"

// #include <core/matrixFreePDE.h>
#include <mpi.h>
#include <vector>

TEST_CASE("Compute invM")
{
  SECTION("1D Subdivided Hyper Rectangle")
  {
    SECTION("Scalar 1st Degree")
    {
      unsigned int dim    = 1;
      unsigned int degree = 1;

      // Initialize triangulation for a unit cube
      std::vector<unsigned int> subdivisions(dim, 1);

      dealii::Triangulation<1> triangulation;
      dealii::GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                        subdivisions,
                                                        dealii::Point<1>(),
                                                        dealii::Point<1>(1.0));

      // Create scalar finite element space
      dealii::FE_Q<1>     fe_q(dealii::QGaussLobatto<1>(degree + 1));
      dealii::FESystem<1> fe_system(fe_q, 1);

      // Create DoFHandler
      dealii::DoFHandler<1> dof_handler(triangulation);
      dof_handler.distribute_dofs(fe_system);

      REQUIRE(dof_handler.n_dofs() == 2);
    }
  }
  SECTION("2D")
  {}
  SECTION("3D")
  {}
}