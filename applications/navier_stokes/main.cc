// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/dependencies.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/problem.h>
#include <prismspf/core/solve_block.h>
#include <prismspf/core/type_enums.h>

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

  // TODO: Add documentation
  std::vector<FieldAttributes> fields = {FieldAttributes("u", Vector),
                                         FieldAttributes("u_star", Vector),
                                         FieldAttributes("p"),
                                         FieldAttributes("p_star"),
                                         FieldAttributes("phi")};

  SolveBlock diffusion;
  diffusion.id            = 2;
  diffusion.field_indices = {0}; // u
  diffusion.solve_type    = Linear;
  diffusion.solve_timing  = Initialized;
  diffusion.dependencies_rhs =
    make_dependency_set(fields, {"old_1(u)", "old_2(u)", "grad(old_1(p_star))"});
  diffusion.dependencies_lhs = make_dependency_set(
    fields,
    {"lhs(u)", "grad(lhs(u))", "old_1(u_star)", "grad(old_1(u_star))"});

  SolveBlock projection;
  projection.id               = 3;
  projection.field_indices    = {4}; // phi
  projection.solve_type       = Linear;
  projection.solve_timing     = Uninitialized;
  projection.dependencies_rhs = make_dependency_set(fields, {"grad(u)"});
  projection.dependencies_lhs = make_dependency_set(fields, {"grad(lhs(phi))"});

  SolveBlock pressure_correction;
  pressure_correction.id               = 4;
  pressure_correction.field_indices    = {2}; // p
  pressure_correction.solve_type       = Explicit;
  pressure_correction.solve_timing     = Initialized;
  pressure_correction.dependencies_rhs = make_dependency_set(fields, {"old_1(p)", "phi"});

  SolveBlock extrapolation;
  extrapolation.id            = 1;
  extrapolation.field_indices = {1, 3}; // u_star & p_star
  extrapolation.solve_type    = Explicit;
  extrapolation.solve_timing  = Uninitialized;
  extrapolation.dependencies_rhs =
    make_dependency_set(fields, {"u", "old_1(u)", "p", "phi", "old_1(phi)"});

  std::vector<SolveBlock> solve_blocks(
    {diffusion, projection, pressure_correction, extrapolation});

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
