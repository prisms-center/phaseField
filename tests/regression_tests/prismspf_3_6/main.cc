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
  prisms::MPI_InitFinalize mpi_init(argc, argv, 1);

  // Restrict deal.II console printing
  dealii::deallog.depth_console(0);

  // Parse the command line options (if there are any) to get the name of the input
  // file
  ParseCMDOptions cli_options(argc, argv);
  std::string     parameters_filename = cli_options.get_parameters_filename();

  constexpr unsigned int dim    = 2;
  constexpr unsigned int degree = 2;

  std::vector<FieldAttributes> field_attributes = {
    FieldAttributes("n", Scalar),
    FieldAttributes("mg_n", Scalar),
    FieldAttributes("f_tot", Scalar),
  };
  std::vector<SolveBlock> solve_blocks;
  SolveBlock              exp_block(1,
                       Explicit,
                       Primary,
                                    {0},
                       make_dependency_set(field_attributes,
                                                        {"old_1(n)", "grad(old_1(n))"}));
  SolveBlock              pp_block(2,
                      Explicit,
                      PostProcess,
                                   {1, 2},
                      make_dependency_set(field_attributes, {"n", "grad(n)"}));
  solve_blocks.push_back(exp_block);
  solve_blocks.push_back(pp_block);

  UserInputParameters<dim>       user_inputs(parameters_filename);
  PhaseFieldTools<dim>           pf_tools;
  CustomPDE<dim, degree, double> pde_operator(user_inputs, pf_tools);
  Problem<dim, degree, double>   problem(field_attributes,
                                       solve_blocks,
                                       user_inputs,
                                       pf_tools,
                                       pde_operator);
  problem.solve();

  return 0;
}
