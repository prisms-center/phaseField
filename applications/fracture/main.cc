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

  constexpr unsigned int dim    = 3;
  constexpr unsigned int degree = 1;

  std::vector<FieldAttributes> fields = {
    FieldAttributes("n", Scalar), // FieldIndex::n     = 0  phase field (explicit)
    FieldAttributes("u", Vector), // FieldIndex::u     = 1  displacement (linear)
    FieldAttributes("dndt",
                    Scalar), // FieldIndex::dndt  = 2  crack driving force (auxiliary)
    FieldAttributes("Ex", Scalar), // FieldIndex::Ex    = 3  stiffness mask (constant)
    FieldAttributes("Gx", Scalar), // FieldIndex::Gx    = 4  toughness mask (constant)
    FieldAttributes("f_tot",
                    Scalar), // FieldIndex::f_tot = 5  total energy density (postprocess)
    FieldAttributes("s11", Scalar), // FieldIndex::s11   = 6  stress (postprocess)
    FieldAttributes("s12", Scalar), // FieldIndex::s12   = 7  stress (postprocess)
    FieldAttributes("s22", Scalar), // FieldIndex::s22   = 8  stress (postprocess)
    FieldAttributes("f_int",
                    Scalar), // FieldIndex::f_int = 9  interfacial energy (postprocess)
    FieldAttributes("f_el",
                    Scalar), // FieldIndex::f_el  = 10 elastic energy (postprocess)
    FieldAttributes("s33",
                    Scalar), // FieldIndex::s33   = 11 stress (postprocess; 0 at dim==2)
    FieldAttributes("s13",
                    Scalar), // FieldIndex::s13   = 12 stress (postprocess; 0 at dim==2)
    FieldAttributes("s23",
                    Scalar), // FieldIndex::s23   = 13 stress (postprocess; 0 at dim==2)
  };

  // Block 0: constant mask fields: initialized once from ICs
  SolveBlock const_block;
  const_block.id            = 0;
  const_block.solve_type    = Constant;
  const_block.solve_timing  = Initialized;
  const_block.field_indices = {idx(FieldIndex::Ex), idx(FieldIndex::Gx)};

  // Block 1: explicit n update using previous-step n and dndt
  SolveBlock n_block;
  n_block.id               = 1;
  n_block.solve_type       = Explicit;
  n_block.solve_timing     = Primary;
  n_block.field_indices    = {idx(FieldIndex::n)};
  n_block.dependencies_rhs = make_dependency_set(fields, {"old_1(n)", "old_1(dndt)"});

  // Block 2: linear u solve: driven by analytical Dirichlet BCs, no body force
  SolveBlock u_block;
  u_block.id               = 2;
  u_block.solve_type       = Linear;
  u_block.solve_timing     = Secondary;
  u_block.field_indices    = {idx(FieldIndex::u)};
  u_block.dependencies_rhs = make_dependency_set(fields, {});
  u_block.dependencies_lhs = make_dependency_set(fields, {"grad(lhs(u))", "n", "Ex"});

  // Block 3: auxiliary dndt: computed each step from current n, u, Ex, Gx
  SolveBlock dndt_block;
  dndt_block.id            = 3;
  dndt_block.solve_type    = Explicit;
  dndt_block.solve_timing  = Secondary;
  dndt_block.field_indices = {idx(FieldIndex::dndt)};
  dndt_block.dependencies_rhs =
    make_dependency_set(fields, {"n", "grad(n)", "grad(u)", "Ex", "Gx"});

  // Block 4: postprocessed fields: computed only on output steps
  SolveBlock pp_block;
  pp_block.id            = 4;
  pp_block.solve_type    = Explicit;
  pp_block.solve_timing  = PostProcess;
  pp_block.field_indices = {
    idx(FieldIndex::f_tot),
    idx(FieldIndex::s11),
    idx(FieldIndex::s12),
    idx(FieldIndex::s22),
    idx(FieldIndex::f_int),
    idx(FieldIndex::f_el),
    idx(FieldIndex::s33),
    idx(FieldIndex::s13),
    idx(FieldIndex::s23),
  };
  pp_block.dependencies_rhs =
    make_dependency_set(fields, {"n", "grad(n)", "grad(u)", "Ex", "Gx"});

  std::vector<SolveBlock> solve_blocks(
    {const_block, n_block, u_block, dndt_block, pp_block});

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
