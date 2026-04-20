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

  /**
   * We have four fields in this application.
   *   U - The dimensionless supersaturation
   *   phi - The solid/liquid order parameter
   *   xi - The auxiliary field used to split the order parameter evolution equation
   *   c - The concentration.
   *
   * The first three equations are explicit with U and phi evolving with a forward Euler
   * time integration scheme. The interesting particle of the equation is xi, which we use
   * to make the evaluation of phi easier. This auxiliary field, xi, has no initial
   * condition and is only used to evolve the order parameter.
   *
   * The last field is the concentration, which we only need for postprocessing.
   */
  std::vector<FieldAttributes> fields = {FieldAttributes("U"),
                                         FieldAttributes("phi"),
                                         FieldAttributes("xi"),
                                         FieldAttributes("c")};

  SolveBlock explicits;
  explicits.id               = 0;
  explicits.solve_type       = Explicit;
  explicits.solve_timing     = Initialized;
  explicits.field_indices    = {0, 1};
  explicits.dependencies_rhs = make_dependency_set(
    fields,
    {"old_1(U)", "grad(old_1(U))", "old_1(phi)", "grad(old_1(phi))", "old_1(xi)"});

  SolveBlock xi_solve;
  xi_solve.id               = 1;
  xi_solve.solve_type       = Explicit;
  xi_solve.solve_timing     = Uninitialized;
  xi_solve.field_indices    = {2};
  xi_solve.dependencies_rhs = make_dependency_set(fields, {"U", "phi", "grad(phi)"});

  SolveBlock pp_solve;
  pp_solve.id               = 2;
  pp_solve.solve_type       = Explicit;
  pp_solve.solve_timing     = PostProcess;
  pp_solve.field_indices    = {3};
  pp_solve.dependencies_rhs = make_dependency_set(fields, {"U", "phi"});

  std::vector<SolveBlock> solves({explicits, xi_solve, pp_solve});

  UserInputParameters<dim>       user_inputs(cli_options.get_parameters_filename());
  PhaseFieldTools<dim>           pf_tools;
  CustomPDE<dim, degree, double> pde_operator(user_inputs, pf_tools);
  Problem<dim, degree, double>   problem(fields,
                                       solves,
                                       user_inputs,
                                       pf_tools,
                                       pde_operator);
  problem.solve();

  return 0;
}
