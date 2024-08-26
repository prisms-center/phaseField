/// Header files
#include "customPDE.h"

#include "ICs_and_BCs.cc"
#include "equations.cc"

#include "../../include/ParseCommandLineOpts.h"
#include "../../include/inputFileReader.h"
#include "../../src/variableAttributeLoader/variableAttributeLoader.cc"

// Header file for postprocessing that may or may not exist
#ifdef POSTPROCESS_FILE_EXISTS
#  include "postprocess.cc"
#else
void
variableAttributeLoader::loadPostProcessorVariableAttributes()
{}
#endif

// Header files for nucleation that may or may not exist
#ifdef NUCLEATION_FILE_EXISTS
#  include "nucleation.cc"

#  include <random>
#  include <time.h>
#endif

// main
int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize
    mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);

  // Parse the command line options (if there are any) to get the name of the
  // input file
  std::string parameters_filename;
  try
    {
      ParseCommandLineOpts cli_options(argc, argv);
      parameters_filename = cli_options.getParametersFilename();
    }
  catch (const char *msg)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------" << std::endl;
      std::cerr << "PRISMS-PF: Exception on processing: " << std::endl
                << msg << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------" << std::endl;
      return 1;
    }

  // Run the main part of the code
  try
    {
      dealii::deallog.depth_console(0);

      // Before fully parsing the parameter file, we need to know how many field
      // variables there are and whether they are scalars or vectors, how many
      // postprocessing variables there are, how many sets of elastic constants
      // there are, and how many user-defined constants there are.

      variableAttributeLoader variable_attributes;
      inputFileReader         input_file_reader(parameters_filename, variable_attributes);

      // Continue based on the number of dimensions and degree of the elements
      // specified in the input file
      switch (input_file_reader.number_of_dimensions)
        {
          case 2:
            {
              userInputParameters<2> userInputs(input_file_reader,
                                                input_file_reader.parameter_handler,
                                                variable_attributes);
              switch (userInputs.degree)
                {
                  case (1):
                    {
                      customPDE<2, 1> problem(userInputs);
                      problem.buildFields();
                      problem.init();
                      problem.solve();
                      break;
                    }
                  case (2):
                    {
                      customPDE<2, 2> problem(userInputs);
                      problem.buildFields();
                      problem.init();
                      problem.solve();
                      break;
                    }
                  case (3):
                    {
                      customPDE<2, 3> problem(userInputs);
                      problem.buildFields();
                      problem.init();
                      problem.solve();
                      break;
                    }
                }
              break;
            }
          case 3:
            {
              userInputParameters<3> userInputs(input_file_reader,
                                                input_file_reader.parameter_handler,
                                                variable_attributes);
              switch (userInputs.degree)
                {
                  case (1):
                    {
                      customPDE<3, 1> problem(userInputs);
                      problem.buildFields();
                      problem.init();
                      problem.solve();
                      break;
                    }
                  case (2):
                    {
                      customPDE<3, 2> problem(userInputs);
                      problem.buildFields();
                      problem.init();
                      problem.solve();
                      break;
                    }
                  case (3):
                    {
                      customPDE<3, 3> problem(userInputs);
                      problem.buildFields();
                      problem.init();
                      problem.solve();
                      break;
                    }
                }
            }
            break;
        }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------" << std::endl;
      std::cerr << "PRISMS-PF: Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------" << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------" << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------" << std::endl;
      return 1;
    }

  return 0;
}
