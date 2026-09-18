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
  constexpr unsigned int number_of_fields = 6;

  std::vector<FieldAttributes> field_attributes;

  for (unsigned int i = 0; i < number_of_fields; i++)
    {
      std::string field_name = "n" + std::to_string(i);
      FieldAttributes attr(field_name, Scalar);
      // To mark grains for remapping, give them a grain_reassignment_block_id
      // Fields with the same grain_reassignment_block_id will be grouped together
      // and treated as a set of order parameters for grain reassignment.
      // Then in the input parameters file, set use grain remapping = true
      attr.grain_reassignment_block_id = 1;
      field_attributes.push_back(attr);
    }

  field_attributes.emplace_back("sum2op", Scalar);
  field_attributes.emplace_back("F", Scalar);
  field_attributes.emplace_back("max_op_id", Scalar);
  
  SolveBlock exp_block;
  exp_block.id            = 1;
  exp_block.solve_type    = Explicit;
  exp_block.solve_timing  = Primary;

  // The order parameters depend on the old value and gradient of the order parameters
  const Dependency old_1_val_and_grad(EvalFlags::nothing, EvalFlags::nothing,
    {EvalFlags::values | EvalFlags::gradients});
  for (unsigned int i = 0; i < number_of_fields; i++)
    {
      exp_block.field_indices.insert(i);
      exp_block.dependencies_rhs[i] = old_1_val_and_grad;
    }

  SolveBlock pp_block;
  pp_block.id               = 2;
  pp_block.solve_type       = Explicit;
  pp_block.solve_timing     = PostProcess;
  pp_block.field_indices    = {number_of_fields, number_of_fields + 1, number_of_fields + 2};
  
  // The postprocessing block depends on the current values and gradients of the order parameters
  const Dependency current_val_and_grad(EvalFlags::values | EvalFlags::gradients);
  for (unsigned int i = 0; i < number_of_fields; i++)
    {
      pp_block.dependencies_rhs[i] = current_val_and_grad;
    }

  std::vector<SolveBlock> solve_blocks({exp_block, pp_block});

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
