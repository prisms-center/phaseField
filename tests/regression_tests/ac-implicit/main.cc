// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/problem.h>

using namespace prisms;

int
main(int argc, char *argv[])
{
  try
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

      constexpr unsigned int dim    = 2;
      constexpr unsigned int degree = 2;

      std::vector<FieldAttributes> field_attributes = {
        FieldAttributes("n", Scalar),
        FieldAttributes("mg_n", Scalar),
        FieldAttributes("f_tot", Scalar),
      };
      std::vector<SolveGroup> solve_groups;
      SolveGroup              newton_group(
        1,
        Newton,
        Primary,
        {0},
        make_dependency_set(field_attributes, {"n", "grad(n)", "old_1(n)"}),
        make_dependency_set(field_attributes, {"n", "change(n)", "grad(change(n))"}));
      SolveGroup pp_group(2,
                          Explicit,
                          PostProcess,
                          {1, 2},
                          make_dependency_set(field_attributes, {"n", "grad(n)"}));
      solve_groups.push_back(newton_group);
      solve_groups.push_back(pp_group);

      UserInputParameters<dim>       user_inputs(parameters_filename);
      PhaseFieldTools<dim>           pf_tools;
      CustomPDE<dim, degree, double> pde_operator(user_inputs, pf_tools);
      Problem<dim, degree, double>   problem(field_attributes,
                                           solve_groups,
                                           user_inputs,
                                           pf_tools,
                                           pde_operator);
      problem.solve();
    }

  catch (std::exception &exc)
    {
      std::cerr << '\n'
                << '\n'
                << "----------------------------------------------------" << '\n';
      std::cerr << "Exception on: " << '\n'
                << exc.what() << '\n'
                << "Aborting!" << '\n'
                << "----------------------------------------------------" << '\n';
      return 1;
    }

  catch (...)
    {
      std::cerr << '\n'
                << '\n'
                << "----------------------------------------------------" << '\n';
      std::cerr << "Unknown exception!" << '\n'
                << "Aborting!" << '\n'
                << "----------------------------------------------------" << '\n';
      return 1;
    }

  return 0;
}
