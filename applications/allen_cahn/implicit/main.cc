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
  constexpr unsigned int degree = 1;

  std::vector<FieldAttributes> field_attributes = {
    FieldAttributes("n", Scalar),
    FieldAttributes("mg_n", Scalar),
    FieldAttributes("f_tot", Scalar),
  };

  SolveBlock newton_block;
  newton_block.id            = 1;
  newton_block.solve_type    = Newton;
  newton_block.solve_timing  = Primary;
  newton_block.field_indices = {0};
  newton_block.dependencies_rhs =
    make_dependency_set(field_attributes, {"n", "grad(n)", "old_1(n)"});
  newton_block.dependencies_lhs =
    make_dependency_set(field_attributes, {"n", "change(n)", "grad(change(n))"});

  SolveBlock pp_block;
  pp_block.id               = 2;
  pp_block.solve_type       = Explicit;
  pp_block.solve_timing     = PostProcess;
  pp_block.field_indices    = {1, 2};
  pp_block.dependencies_rhs = make_dependency_set(field_attributes, {"n", "grad(n)"});

  std::vector<SolveBlock> solve_blocks({newton_block, pp_block});

  UserInputParameters<dim>       user_inputs(cli_options.get_parameters_filename());
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