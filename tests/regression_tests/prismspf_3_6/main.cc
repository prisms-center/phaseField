// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/field_attributes.h>
#include <prismspf/core/parse_cmd_options.h>
#include <prismspf/core/phase_field_tools.h>
#include <prismspf/core/problem.h>
#include <prismspf/core/solve_group.h>

#include <prismspf/user_inputs/input_file_reader.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

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
      prisms::ParseCMDOptions cli_options(argc, argv);
      std::string             parameters_filename = cli_options.get_parameters_filename();

      constexpr unsigned int dim    = 1;
      constexpr unsigned int degree = 1;

      std::vector<prisms::FieldAttributes> field_attributes = {
        prisms::FieldAttributes("n", prisms::FieldInfo::TensorRank::Scalar),
        prisms::FieldAttributes("mg_n", prisms::FieldInfo::TensorRank::Scalar),
        prisms::FieldAttributes("f_tot", prisms::FieldInfo::TensorRank::Scalar),
      };
      std::vector<prisms::SolveGroup> solve_groups;

      prisms::UserInputParameters<dim>                              user_inputs();
      prisms::PhaseFieldTools<dim>                                  pf_tools;
      std::shared_ptr<prisms::PDEOperatorBase<dim, degree, double>> pde_operator =
        std::make_shared<prisms::CustomPDE<dim, degree, double>>(user_inputs, pf_tools);
      prisms::Problem<dim, degree, double> problem(user_inputs, pf_tools, pde_operator);
      problem.run();
    }

  catch (std::exception &exc)
    {
      std::cerr << '\n'
                << '\n'
                << "----------------------------------------------------" << '\n';
      std::cerr << "Exception on processing: " << '\n'
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
