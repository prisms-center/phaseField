// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/problem.h>

using namespace prismspf;

int
main(int argc, char *argv[])
{
  prismspf::MPIInitFinalize mpi_init(argc, argv);
  ParseCMDOptions           cli_options(argc, argv);

  constexpr unsigned int dim    = 2;
  constexpr unsigned int degree = 1;

  std::vector<FieldAttributes> fields = {
    FieldAttributes("n",     Scalar), // 0 - phase field to track fracture (scalar, explicit)
    FieldAttributes("u",     Vector), // 1 - displacement (vector, linear)
    FieldAttributes("dndt",  Scalar), // 2 - crack driving force (scalar, auxiliary)
    FieldAttributes("Ex",    Scalar), // 3 - stiffness mask (scalar, constant)
    FieldAttributes("Gx",    Scalar), // 4 - toughness mask (scalar, constant)
    FieldAttributes("f_tot", Scalar), // 5 - total energy density (scalar, postprocess)
    FieldAttributes("s11",   Scalar), // 6 - stress component (scalar, postprocess)
    FieldAttributes("s12",   Scalar), // 7 - stress component (scalar, postprocess)
    FieldAttributes("s22",   Scalar), // 8 - stress component (scalar, postprocess)
    FieldAttributes("f_int", Scalar), // 9 - interfacial energy density (scalar, postprocess)
    FieldAttributes("f_el",  Scalar), // 10 - elastic energy density (scalar, postprocess)
  };

  // Block 0: explicit n update using previous-step n and dndt
  SolveBlock n_block;
  n_block.id            = 0;
  n_block.solve_type    = Explicit;
  n_block.solve_timing  = Primary;
  n_block.field_indices = {0};
  n_block.dependencies_rhs =
    make_dependency_set(fields, {"old_1(n)", "old_1(dndt)"});

  // Block 1: constant mask fields — initialized once from ICs
  SolveBlock const_block;
  const_block.id            = 1;
  const_block.solve_type    = Constant;
  const_block.solve_timing  = Initialized;
  const_block.field_indices = {3, 4};

  // Block 2: linear u solve — driven by analytical Dirichlet BCs, no body force
  SolveBlock u_block;
  u_block.id            = 2;
  u_block.solve_type    = Linear;
  u_block.solve_timing  = Secondary;
  u_block.field_indices = {1};
  u_block.dependencies_rhs = make_dependency_set(fields, {});
  u_block.dependencies_lhs =
    make_dependency_set(fields, {"grad(lhs(u))", "n", "Ex"});

  // Block 3: auxiliary dndt — computed each step from current n, u, Ex, Gx
  SolveBlock dndt_block;
  dndt_block.id            = 3;
  dndt_block.solve_type    = Explicit;
  dndt_block.solve_timing  = Secondary;
  dndt_block.field_indices = {2};
  dndt_block.dependencies_rhs =
    make_dependency_set(fields, {"n", "grad(n)", "grad(u)", "Ex", "Gx"});

  // Block 4: postprocessed fields — computed only on output steps
  SolveBlock pp_block;
  pp_block.id            = 4;
  pp_block.solve_type    = Explicit;
  pp_block.solve_timing  = PostProcess;
  pp_block.field_indices = {5, 6, 7, 8, 9, 10};
  pp_block.dependencies_rhs =
    make_dependency_set(fields, {"n", "grad(n)", "grad(u)", "Ex", "Gx"});

  std::vector<SolveBlock> solve_blocks(
    {n_block, const_block, u_block, dndt_block, pp_block});

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
