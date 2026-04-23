// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/problem.h>

using namespace prisms;

int
main(int argc, char *argv[])
{
  // Initialize MPI
  prisms::MPIInitFinalize mpi_init(argc, argv);

  // Parse the command line options (if there are any) to get the name of the input
  // file
  ParseCMDOptions cli_options(argc, argv);

  constexpr unsigned int dim    = 2;
  constexpr unsigned int degree = 2;

  std::vector<FieldAttributes> fields = {
    FieldAttributes("c"),
    FieldAttributes("n"),
    FieldAttributes("nucleation_rate", Scalar, true, {1}),
  };

  SolveBlock explicits;
  explicits.id            = 0;
  explicits.solve_type    = Explicit;
  explicits.solve_timing  = Primary;
  explicits.field_indices = {0, 1};
  explicits.dependencies_rhs =
    make_dependency_set(fields,
                        {"old_1(c)", "grad(old_1(c))", "old_1(n)", "grad(old_1(n))"});

  SolveBlock nucleation;
  nucleation.id               = 1;
  nucleation.solve_type       = Explicit;
  nucleation.solve_timing     = NucleationRate;
  nucleation.field_indices    = {2};
  nucleation.dependencies_rhs = make_dependency_set(fields, {"n", "c"});

  std::vector<SolveBlock> solve_blocks({explicits, nucleation});

  UserInputParameters<dim>       user_inputs(cli_options.get_parameters_filename());
  PhaseFieldTools<dim>           pf_tools;
  CustomPDE<dim, degree, double> pde_operator(user_inputs, pf_tools);
  Problem<dim, degree, double>   problem(fields,
                                       solve_blocks,
                                       user_inputs,
                                       pf_tools,
                                       pde_operator);
  problem.solve();

  return 0;
}
