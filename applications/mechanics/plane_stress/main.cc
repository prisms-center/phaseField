// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/problem.h>

using namespace prismspf;

int
main(int argc, char *argv[])
{
  // Initialize MPI
  prismspf::MPIInitFinalize mpi_init(argc, argv);

  // Parse the command line options (if there are any) to get the name of the input
  // file
  ParseCMDOptions cli_options(argc, argv);

  constexpr unsigned int dim    = 2;
  constexpr unsigned int degree = 2;

  std::vector<FieldAttributes> fields = {
    FieldAttributes("u", Vector),
    FieldAttributes("epsilon_xx", Scalar),
    FieldAttributes("epsilon_yy", Scalar),
    FieldAttributes("gamma_xy", Scalar),
    FieldAttributes("sigma_xx", Scalar),
    FieldAttributes("sigma_yy", Scalar),
    FieldAttributes("sigma_xy", Scalar),
  };

  SolveBlock linear_solve;
  linear_solve.id               = 1;
  linear_solve.solve_type       = Linear;
  linear_solve.solve_timing     = Uninitialized;
  linear_solve.field_indices    = {0};
  linear_solve.dependencies_lhs = make_dependency_set(fields, {"grad(lhs(u))"});

  SolveBlock pp_block;
  pp_block.id         = 2;
  pp_block.solve_type = Explicit, pp_block.solve_timing = PostProcess,
  pp_block.field_indices    = {1, 2, 3, 4, 5, 6},
  pp_block.dependencies_rhs = make_dependency_set(fields, {"grad(u)"});

  std::vector<SolveBlock> solve_blocks({linear_solve, pp_block});

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
