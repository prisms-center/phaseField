#include "custom_pde.h"

#include <core/parse_cmd_options.h>
#include <core/pde_problem.h>
#include <core/user_inputs/input_file_reader.h>
#include <core/user_inputs/user_input_parameters.h>
#include <core/variable_attribute_loader.h>
#include <core/variable_attributes.h>

int
main(int argc, char *argv[])
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);

      // Parse the command line options (if there are any) to get the name of the input
      // file
      parseCMDOptions cli_options(argc, argv);
      std::string     parameters_filename = cli_options.get_parameters_filename();

      // Restrict deal.II console printing
      dealii::deallog.depth_console(0);

      // Before fully parsing the parameter file, we need to know how many field
      // variables there are and whether they are scalars or vectors, how many
      // postprocessing variables there are, how many sets of elastic constants
      // there are, and how many user-defined constants there are.
      //
      // This is done with the derived class of `variableAttributeLoader`,
      // `customAttributeLoader`.
      customAttributeLoader attribute_loader;
      attribute_loader.init_variable_attributes();
      AttributesList var_attributes = attribute_loader.get_var_attributes();
      AttributesList pp_attributes  = attribute_loader.get_pp_attributes();

      // Load in parameters
      inputFileReader input_file_reader(parameters_filename,
                                        var_attributes,
                                        pp_attributes);

      // Run problem based on the number of dimensions and element degree
      switch (input_file_reader.number_of_dimensions)
        {
          case 1:
            {
              userInputParameters<1> user_inputs(input_file_reader,
                                                 input_file_reader.parameter_handler);
              switch (user_inputs.spatial_discretization.degree)
                {
                  case 1:
                    {
                      PDEProblem<1, 1> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 2:
                    {
                      PDEProblem<1, 2> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 3:
                    {
                      PDEProblem<1, 3> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 4:
                    {
                      PDEProblem<1, 4> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 5:
                    {
                      PDEProblem<1, 5> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 6:
                    {
                      PDEProblem<1, 6> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  default:
                    throw std::runtime_error("Invalid element degree");
                }
              break;
            }
          case 2:
            {
              userInputParameters<2> user_inputs(input_file_reader,
                                                 input_file_reader.parameter_handler);
              switch (user_inputs.spatial_discretization.degree)
                {
                  case 1:
                    {
                      PDEProblem<2, 1> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 2:
                    {
                      PDEProblem<2, 2> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 3:
                    {
                      PDEProblem<2, 3> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 4:
                    {
                      PDEProblem<2, 4> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 5:
                    {
                      PDEProblem<2, 5> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 6:
                    {
                      PDEProblem<2, 6> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  default:
                    throw std::runtime_error("Invalid element degree");
                }
              break;
            }
          case 3:
            {
              userInputParameters<3> user_inputs(input_file_reader,
                                                 input_file_reader.parameter_handler);
              switch (user_inputs.spatial_discretization.degree)
                {
                  case 1:
                    {
                      PDEProblem<3, 1> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 2:
                    {
                      PDEProblem<3, 2> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 3:
                    {
                      PDEProblem<3, 3> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 4:
                    {
                      PDEProblem<3, 4> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 5:
                    {
                      PDEProblem<3, 5> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  case 6:
                    {
                      PDEProblem<3, 6> problem(user_inputs);
                      problem.run();
                      break;
                    }
                  default:
                    throw std::runtime_error("Invalid element degree");
                }
              break;
            }
          default:
            throw std::runtime_error("Invalid number of dimensions");
        }
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