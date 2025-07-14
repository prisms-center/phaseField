// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/user_inputs/load_initial_condition_parameters.h>

#include <prismspf/config.h>

#include "catch.hpp"

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Test the initial condition read in.
 */
TEST_CASE("Initial condition read in")
{
  SECTION("Initial condition read in disallowed")
  {
    InitialConditionFile           file;
    LoadInitialConditionParameters parameters;
    parameters.add_initial_condition_file(file);
    REQUIRE(parameters.get_n_initial_condition_files() == 0);
    parameters.clear();
  }
  SECTION("Mismatching number of file and simulation variable names")
  {
    InitialConditionFile file;
    file.file_variable_names       = {"n", "n1"};
    file.simulation_variable_names = {"n"};

    LoadInitialConditionParameters parameters;
    parameters.set_read_initial_conditions_from_file(true);
    parameters.add_initial_condition_file(file);
    REQUIRE_THROWS(parameters.postprocess_and_validate());
    parameters.clear();
  }
}

PRISMS_PF_END_NAMESPACE