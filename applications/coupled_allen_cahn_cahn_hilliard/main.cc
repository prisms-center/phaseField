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
  dealii::Utilities::MPI::MPI_InitFinalize
    mpi_init(argc, argv, dealii::numbers::invalid_unsigned_int);

  // Restrict deal.II console printing
  dealii::deallog.depth_console(0);

  // Parse the command line options (if there are any) to get the name of the input
  // file
  ParseCMDOptions cli_options(argc, argv);

  constexpr unsigned int dim    = 2;
  constexpr unsigned int degree = 1;

  std::vector<FieldAttributes> fields = {
    FieldAttributes("c", Scalar),
    FieldAttributes("n", Scalar),
    FieldAttributes("f_tot", Scalar),
    FieldAttributes("mag_grad_c", Scalar),

  };

  SolveBlock exp_block(
    0,
    Explicit,
    Primary,
    {0, 1},
    make_dependency_set(fields,
                        {"old_1(c)", "grad(old_1(c))", "old_1(n)", "grad(old_1(n))"}));
  SolveBlock              pp_block(1,
                      Explicit,
                      PostProcess,
                                   {2, 3},
                      make_dependency_set(fields, {"n", "grad(n)", "c", "grad(c)"}));
  std::vector<SolveBlock> solve_blocks({exp_block, pp_block});

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
