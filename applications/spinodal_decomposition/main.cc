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
  std::string     parameters_filename = cli_options.get_parameters_filename();

  constexpr unsigned int dim    = 2; // TODO change to 3 (original app)
  constexpr unsigned int degree = 2; // TODO change to 1 (original app)

  std::vector<FieldAttributes> fields = {FieldAttributes("c"),
                                         FieldAttributes("mu"),
                                         FieldAttributes("f_tot")};
  std::vector<SolveGroup>      solve_groups;
  SolveGroup                   c_group;
  c_group.id               = 0;
  c_group.solve_type       = Explicit;
  c_group.solve_timing     = Initialized;
  c_group.field_indices    = {0};
  c_group.dependencies_rhs = make_dependency_set(fields, {"old_1(c)", "grad(old_1(mu))"});

  SolveGroup mu_group;
  mu_group.id               = 1;
  mu_group.solve_type       = Explicit;
  mu_group.solve_timing     = Uninitialized;
  mu_group.field_indices    = {1};
  mu_group.dependencies_rhs = make_dependency_set(fields, {"c", "grad(c)"});

  SolveGroup pp_group;
  pp_group.id               = 2;
  pp_group.solve_type       = Explicit;
  pp_group.solve_timing     = PostProcess;
  pp_group.field_indices    = {2};
  pp_group.dependencies_rhs = make_dependency_set(fields, {"c", "grad(c)"});

  solve_groups.push_back(c_group);
  solve_groups.push_back(mu_group);
  solve_groups.push_back(pp_group);

  UserInputParameters<dim>       user_inputs(parameters_filename);
  PhaseFieldTools<dim>           pf_tools;
  CustomPDE<dim, degree, double> pde_operator(user_inputs, pf_tools);
  Problem<dim, degree, double>   problem(fields,
                                       solve_groups,
                                       user_inputs,
                                       pf_tools,
                                       pde_operator);
  problem.solve();

  return 0;
}
