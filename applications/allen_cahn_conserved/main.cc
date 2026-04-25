// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/problem.h>
#include <prismspf/core/solve_block.h>
#include <prismspf/core/type_enums.h>

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
    FieldAttributes("c", Scalar),
    FieldAttributes("mu", Scalar),
    FieldAttributes("dummy", Scalar),
    FieldAttributes("mg_n", Scalar),
    FieldAttributes("f_tot", Scalar),
  };

  SolveBlock c_block;
  c_block.id            = 0;
  c_block.solve_type    = Explicit;
  c_block.solve_timing  = Primary;
  c_block.field_indices = {0};
  c_block.dependencies_rhs =
    make_dependency_set(field_attributes, {"old_1(c)", "old_1(mu)"});

  SolveBlock mu_block;
  mu_block.id               = 1;
  mu_block.solve_type       = Explicit;
  mu_block.solve_timing     = Uninitialized;
  mu_block.field_indices    = {1};
  mu_block.dependencies_rhs = make_dependency_set(field_attributes, {"c", "grad(c)"});

  SolveBlock cs_block;
  cs_block.id               = 2;
  cs_block.solve_type       = CustomSolver;
  cs_block.solve_timing     = Uninitialized;
  cs_block.field_indices    = {2};
  cs_block.dependencies_rhs = make_dependency_set(field_attributes, {"mu"});

  SolveBlock pp_block;
  pp_block.id               = 3;
  pp_block.solve_type       = Explicit;
  pp_block.solve_timing     = PostProcess;
  pp_block.field_indices    = {3, 4};
  pp_block.dependencies_rhs = make_dependency_set(field_attributes, {"c", "grad(c)"});

  UserInputParameters<dim>       user_inputs(cli_options.get_parameters_filename());
  PhaseFieldTools<dim>           pf_tools;
  CustomPDE<dim, degree, double> pde_operator(user_inputs, pf_tools);

  // setup CustomSolver
  using CustomFactory = std::function<std::shared_ptr<SolverBase<dim, degree, double>>(
    const SolveBlock &,
    const SolveContext<dim, degree, double> &)>;

  cs_block.custom_solver_factory = CustomFactory(
    [&pde_operator](const SolveBlock &b, const SolveContext<dim, degree, double> &c)
    {
      return std::make_shared<MyCustomSolver<dim, degree, double>>(b, c, pde_operator);
    });

  std::vector<SolveBlock> solve_blocks({c_block, mu_block, cs_block, pp_block});

  Problem<dim, degree, double> problem(field_attributes,
                                       solve_blocks,
                                       user_inputs,
                                       pf_tools,
                                       pde_operator);
  problem.solve();

  return 0;
}
