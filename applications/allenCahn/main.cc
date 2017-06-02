// Allen-Cahn example application

// Header files
#include "parameters.h"
#include "../../include/dealIIheaders.h"
#include "../../include/model_variables.h"
#include "../../include/varBCs.h"
#include "../../include/initialConditions.h"
#include "../../include/matrixFreePDE.h"
#include "customPDE.h"
#include "equations.h"
#include "ICs_and_BCs.h"
#include "postprocess.h"
#include "../../include/initialCondition_template_instantiations.h"
#include "../../include/inputFileReader.h"
#include "../../include/userInputParameters.h"
#include "../../src/userInputParameters/loadUserInputs.cc" // Needs to be included because it contains needs access to the define macros in the preceding files//#include "../../src/userInputParameters/declareInputParameters.cc" // Needs to be included because it contains needs access to the define macros in the preceding files
#include "../../src/userInputParameters/loadInputParameters.cc" // Needs to be included because it contains needs access to the define macros in the preceding files


//main
int main (int argc, char **argv)
{
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,dealii::numbers::invalid_unsigned_int);
    try
    {
        dealii::deallog.depth_console(0);

        dealii::ParameterHandler parameter_handler;
        inputFileReader input_file_reader(parameter_handler,"parameters.in");

        const unsigned int dimension = parameter_handler.get_integer("Number of dimensions");

        switch (dimension)
        {
            case 2:
            {
                userInputParameters<2> userInputs;
                userInputs.loadUserInput();
                userInputs.loadInputParameters(parameter_handler);

                switch (userInputs.degree)
                {
                    case(1):
                    {
                        customPDE<2,1> problem(userInputs);
                        problem.setBCs();
                        problem.buildFields();
                        problem.init ();
                        problem.solve();
                        break;
                    }
                    case(2):
                    {
                        customPDE<2,2> problem(userInputs);
                        problem.setBCs();
                        problem.buildFields();
                        problem.init ();
                        problem.solve();
                        break;
                    }
                    case(3):
                    {
                        customPDE<2,3> problem(userInputs);
                        problem.setBCs();
                        problem.buildFields();
                        problem.init ();
                        problem.solve();
                        break;
                    }
                }
                break;
            }
            case 3:
            {
                userInputParameters<3> userInputs;
                userInputs.loadUserInput();
                userInputs.loadInputParameters(parameter_handler);

                switch (userInputs.degree)
                {
                    case(1):
                    {
                        customPDE<3,1> problem(userInputs);
                        problem.setBCs();
                        problem.buildFields();
                        problem.init ();
                        problem.solve();
                        break;
                    }
                    case(2):
                    {
                        customPDE<3,2> problem(userInputs);
                        problem.setBCs();
                        problem.buildFields();
                        problem.init ();
                        problem.solve();
                        break;
                    }
                    case(3):
                    {
                        customPDE<3,3> problem(userInputs);
                        problem.setBCs();
                        problem.buildFields();
                        problem.init ();
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
        std::cerr << std::endl << std::endl
        << "----------------------------------------------------"
        << std::endl;
        std::cerr << "Exception on processing: " << std::endl
        << exc.what() << std::endl
        << "Aborting!" << std::endl
        << "----------------------------------------------------"
        << std::endl;
        return 1;
    }
    catch (...)
    {
        std::cerr << std::endl << std::endl
        << "----------------------------------------------------"
        << std::endl;
        std::cerr << "Unknown exception!" << std::endl
        << "Aborting!" << std::endl
        << "----------------------------------------------------"
        << std::endl;
        return 1;
    }

    return 0;
}
