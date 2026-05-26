// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/problem.h>

#include <prismspf/user_inputs/constraint_parameters.h>

using namespace prismspf;

int
main(int argc, char *argv[])
{
  // Initialize MPI
  prismspf::MPIInitFinalize mpi_init(argc, argv);

  // Parse the command line options (if there are any) to get the name of the input
  // file
  ParseCMDOptions        cli_options(argc, argv);
  std::string            parameters_filename = cli_options.get_parameters_filename();
  constexpr unsigned int degree              = 1;

  std::vector<FieldAttributes> fields = {FieldAttributes("c"),
                                         FieldAttributes("mu"),
                                         FieldAttributes("f_tot"),
                                         FieldAttributes("Phi")};

  SolveBlock c_block;
  c_block.id               = 0;
  c_block.solve_type       = Explicit;
  c_block.solve_timing     = Initialized;
  c_block.field_indices    = {0};
  c_block.dependencies_rhs = make_dependency_set(fields, {"old_1(c)", "grad(old_1(mu))"});

  SolveBlock mu_block;
  mu_block.id               = 1;
  mu_block.solve_type       = Explicit;
  mu_block.solve_timing     = Uninitialized;
  mu_block.field_indices    = {1};
  mu_block.dependencies_rhs = make_dependency_set(fields, {"c", "grad(c)"});

  SolveBlock Phi_block;
  Phi_block.id               = 3;
  Phi_block.solve_type       = Linear;
  Phi_block.solve_timing     = Uninitialized;
  Phi_block.field_indices    = {3};
  Phi_block.dependencies_rhs = make_dependency_set(fields, {"c"});
  Phi_block.dependencies_lhs = make_dependency_set(fields, {"grad(lhs(Phi))"});

  SolveBlock pp_block;
  pp_block.id               = 2;
  pp_block.solve_type       = Explicit;
  pp_block.solve_timing     = PostProcess;
  pp_block.field_indices    = {2};
  pp_block.dependencies_rhs = make_dependency_set(fields, {"c", "grad(c)", "Phi"});

  std::vector<SolveBlock> solve_blocks({c_block, mu_block, Phi_block, pp_block});

  UserInputParameters<2> user_inputs(parameters_filename);
  user_inputs.spatial_discretization.rectangular_mesh.size =
    dealii::Tensor<1, 2>({100.0, 100.0});
  FieldConstraints<2> &Phi_bcs =
    user_inputs.boundary_parameters.boundary_condition_list["Phi"];
  Phi_bcs.component_constraints[0].conditions[0] = Dirichlet;
  Phi_bcs.component_constraints[0].conditions[1] = Dirichlet;

  // Output at times specified by the benchmark spec.
  for (double output_time : {400})
    {
      auto output_step =
        static_cast<unsigned int>(output_time / user_inputs.temporal_discretization.dt);
      user_inputs.output_parameters.output_list.insert(output_step);
    }
  PhaseFieldTools<2>         pf_tools;
  CustomPDE<degree, double>  pde_operator(user_inputs, pf_tools);
  Problem<2, degree, double> problem(fields,
                                     solve_blocks,
                                     user_inputs,
                                     pf_tools,
                                     pde_operator);
  problem.solve();

  return 0;
}
