// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/user_inputs/nonlinear_solve_parameters.h>

#include <prismspf/config.h>

#include "catch.hpp"

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Test the nonlinear solve parameters.
 */
TEST_CASE("Nonlinear solve parameters")
{
  NonlinearSolveParameters parameters;
  SECTION("Step length")
  {
    NonlinearSolverParameters solver_parameters;

    solver_parameters.step_length = 0.5;
    parameters.set_nonlinear_solve_parameters(0, solver_parameters);
    REQUIRE_NOTHROW(parameters.postprocess_and_validate());
    parameters.clear();

    solver_parameters.step_length = -1.0;
    parameters.set_nonlinear_solve_parameters(0, solver_parameters);
    REQUIRE_THROWS(parameters.postprocess_and_validate());
    parameters.clear();

    solver_parameters.step_length = 0.0;
    parameters.set_nonlinear_solve_parameters(0, solver_parameters);
    REQUIRE_THROWS(parameters.postprocess_and_validate());
    parameters.clear();

    solver_parameters.step_length = 2.0;
    parameters.set_nonlinear_solve_parameters(0, solver_parameters);
    REQUIRE_THROWS(parameters.postprocess_and_validate());
    parameters.clear();
  }
  SECTION("Tolerance value")
  {
    NonlinearSolverParameters solver_parameters;
    solver_parameters.tolerance_value = 1e-10;
    parameters.set_nonlinear_solve_parameters(0, solver_parameters);
    REQUIRE_NOTHROW(parameters.postprocess_and_validate());
    parameters.clear();

    solver_parameters.tolerance_value = -1.0;
    parameters.set_nonlinear_solve_parameters(0, solver_parameters);
    REQUIRE_THROWS(parameters.postprocess_and_validate());
    parameters.clear();

    solver_parameters.tolerance_value = 0.0;
    parameters.set_nonlinear_solve_parameters(0, solver_parameters);
    REQUIRE_THROWS(parameters.postprocess_and_validate());
    parameters.clear();
  }
}

PRISMS_PF_END_NAMESPACE