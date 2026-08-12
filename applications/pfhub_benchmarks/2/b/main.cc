// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
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
  ParseCMDOptions        cli_options(argc, argv);
  std::string            parameters_filename = cli_options.get_parameters_filename();
  constexpr unsigned int degree              = 2;
  constexpr unsigned int p                   = CustomPDE<degree, double>::p;

  std::vector<FieldAttributes> fields;
  fields.reserve(3 + p);
  for (unsigned int i = 0; i < p; ++i)
    {
      FieldAttributes c_attr;
      c_attr.name = "eta_" + std::to_string(i);
      fields.push_back(c_attr);
    }
  fields.push_back(FieldAttributes("c"));
  fields.push_back(FieldAttributes("mu"));
  fields.push_back(FieldAttributes("F"));

  const Dependency old_1_value_and_gradient(EvalFlags::nothing,
                                            EvalFlags::nothing,
                                            {EvalFlags::values | EvalFlags::gradients});
  const Dependency old_1_gradient(EvalFlags::nothing,
                                  EvalFlags::nothing,
                                  {EvalFlags::gradients});
  const Dependency current_value_and_gradient(EvalFlags::values | EvalFlags::gradients);
  const Dependency current_value(EvalFlags::values);

  SolveBlock explicit_block;
  explicit_block.id           = 0;
  explicit_block.solve_type   = Explicit;
  explicit_block.solve_timing = Initialized;
  for (unsigned int i = 0; i < p + 1; ++i)
    {
      explicit_block.field_indices.insert(i);
    }
  for (unsigned int i = 0; i < p; ++i)
    {
      explicit_block.dependencies_rhs[i] = old_1_value_and_gradient;
    }
  explicit_block.dependencies_rhs[p] =
    Dependency(EvalFlags::nothing, EvalFlags::nothing, {EvalFlags::values});
  explicit_block.dependencies_rhs[p + 1] = old_1_gradient;

  SolveBlock mu_block;
  mu_block.id            = 1;
  mu_block.solve_type    = Explicit;
  mu_block.solve_timing  = Uninitialized;
  mu_block.field_indices = {p + 1};
  for (unsigned int i = 0; i < p; ++i)
    {
      mu_block.dependencies_rhs[i] = current_value;
    }
  mu_block.dependencies_rhs[p] = current_value_and_gradient;

  SolveBlock pp_block;
  pp_block.id            = 2;
  pp_block.solve_type    = Explicit;
  pp_block.solve_timing  = PostProcess;
  pp_block.field_indices = {p + 2};
  for (unsigned int i = 0; i < p + 1; ++i)
    {
      pp_block.dependencies_rhs[i] = current_value_and_gradient;
    }

  std::vector<SolveBlock> solve_blocks({explicit_block, mu_block, pp_block});

  UserInputParameters<2> user_inputs(parameters_filename);
  user_inputs.spatial_discretization.rectangular_mesh.size =
    dealii::Tensor<1, 2>({200.0, 200.0});
  // Output at times specified by the benchmark spec.
  for (double output_time : {1000, 10000, 100000, 1000000})
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
