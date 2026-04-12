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
  prisms::MPI_InitFinalize mpi_init(argc, argv, dealii::numbers::invalid_unsigned_int);

  // Restrict deal.II console printing
  dealii::deallog.depth_console(0);

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

  SolveBlock explicits(
    0,
    Explicit,
    Initialized,
    {0, 1},
    make_dependency_set(
      fields,
      {"old_1(U)", "grad(old_1(U))", "old_1(phi)", "grad(old_1(phi))", "old_1(xi)"}));
  SolveBlock xi_solve(1,
                      Explicit,
                      Uninitialized,
                      {2},
                      make_dependency_set(fields, {"U", "phi", "grad(phi)"}));

  SolveBlock pp_solve(2,
                      Explicit,
                      PostProcess,
                      {3},
                      make_dependency_set(fields, {"U", "phi"}));

  std::vector<SolveBlock>        solves({explicits, xi_solve, pp_solve});
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
