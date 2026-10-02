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
    FieldAttributes("n"),         // FieldIndex::n = 0  phase field (explicit)
    FieldAttributes("u", Vector), // FieldIndex::u = 1  displacement (linear)
    FieldAttributes("dndt"),  // FieldIndex::dndt  = 2  crack driving force (auxiliary)
    FieldAttributes("Ex"),    // FieldIndex::Ex    = 3  stiffness mask (constant)
    FieldAttributes("Gx"),    // FieldIndex::Gx    = 4  toughness mask (constant)
    FieldAttributes("f_tot"), // FieldIndex::f_tot = 5  total energy density (postprocess)
    FieldAttributes("f_int"), // FieldIndex::f_int = 6  interfacial energy (postprocess)
    FieldAttributes("f_el"),  // FieldIndex::f_el  = 7 elastic energy (postprocess)
    FieldAttributes("s11"),   // FieldIndex::s11   = 8  stress (postprocess)
  };
  if constexpr (dim > 1)
    {
      fields.emplace_back(FieldAttributes("s12")); // FieldIndex::s33   = 9 stress
      fields.emplace_back(FieldAttributes("s22")); // FieldIndex::s13   = 10 stress
    }
  if constexpr (dim > 2)
    {
      fields.emplace_back(FieldAttributes("s33")); // FieldIndex::s33   = 11 stress
      fields.emplace_back(FieldAttributes("s13")); // FieldIndex::s13   = 12 stress
      fields.emplace_back(FieldAttributes("s23")); // FieldIndex::s23   = 13 stress
    }

  // Block 0: constant mask fields: initialized once from ICs
  SolveBlock const_block;
  const_block.id            = 0;
  const_block.solve_type    = Constant;
  const_block.solve_timing  = Initialized;
  const_block.field_indices = {FieldIndex::Ex, FieldIndex::Gx};

  // Block 1: explicit n update using previous-step n and dndt
  SolveBlock n_block;
  n_block.id               = 1;
  n_block.solve_type       = Explicit;
  n_block.solve_timing     = Primary;
  n_block.field_indices    = {FieldIndex::n};
  n_block.dependencies_rhs = make_dependency_set(fields, {"old_1(n)", "old_1(dndt)"});

  // Block 2: linear u solve: driven by analytical Dirichlet BCs, no body force
  SolveBlock u_block;
  u_block.id               = 2;
  u_block.solve_type       = Linear;
  u_block.solve_timing     = Secondary;
  u_block.field_indices    = {FieldIndex::u};
  u_block.dependencies_rhs = make_dependency_set(fields, {});
  u_block.dependencies_lhs = make_dependency_set(fields, {"grad(lhs(u))", "n", "Ex"});

  // Block 3: auxiliary dndt: computed each step from current n, u, Ex, Gx
  SolveBlock dndt_block;
  dndt_block.id            = 3;
  dndt_block.solve_type    = Explicit;
  dndt_block.solve_timing  = Secondary;
  dndt_block.field_indices = {FieldIndex::dndt};
  dndt_block.dependencies_rhs =
    make_dependency_set(fields, {"n", "grad(n)", "grad(u)", "Ex", "Gx"});

  // Block 4: postprocessed fields: computed only on output steps
  SolveBlock pp_block;
  pp_block.id            = 4;
  pp_block.solve_type    = Explicit;
  pp_block.solve_timing  = PostProcess;
  pp_block.field_indices = {FieldIndex::f_tot,
                            FieldIndex::f_int,
                            FieldIndex::f_el,
                            FieldIndex::s11};
  if constexpr (dim > 1)
    {
      pp_block.field_indices.insert({FieldIndex::s12, FieldIndex::s22});
    }
  if constexpr (dim > 2)
    {
      pp_block.field_indices.insert({FieldIndex::s33, FieldIndex::s13, FieldIndex::s23});
    }

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
