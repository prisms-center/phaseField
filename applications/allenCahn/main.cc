// Allen-Cahn example application

// Header files
#include "../../include/ParseCommandLineOpts.h"
#include "../../include/initialConditions.h"
#include "../../include/initialCondition_template_instantiations.h"
#include "../../include/nonUniformDirichletBC.h"
#include "../../include/nonUniformDirichletBC_template_instantiations.h"
#include "../../include/matrixFreePDE.h"
#include "../../src/models/mechanics/computeStress.h"
#include "../../include/inputFileReader.h"
#include "customPDE.h"
#include "equations.h"
#include "ICs_and_BCs.h"
#include "../../src/variableAttributeLoader/variableAttributeLoader.cc"

// Header file for postprocessing that may or may not exist
#ifdef POSTPROCESS_FILE_EXISTS
#include "postprocess.h"
#else
void variableAttributeLoader::loadPostProcessorVariableAttributes(){}
#endif

// Header files for nucleation that may or may not exist
#ifdef NUCLEATION_FILE_EXISTS
#include <random>
#include <time.h>
#include "nucleation.h"
#endif

//main
int main (int argc, char **argv)
{
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,dealii::numbers::invalid_unsigned_int);

    // Parse the command line options (if there are any) to get the name of the input file
    std::string parameters_filename;
    try
    {
        ParseCommandLineOpts cli_options(argc, argv);
        if (argc == 3){
            if (cli_options.cmdOptionExists("-i")){
                parameters_filename = cli_options.getCmdOption("-i");
                std::cout << "Using the input parameter file: " << parameters_filename << std::endl;
            }
            else {
                throw(0);
            }
        }
        else if (argc == 2){
            if (cli_options.cmdOptionExists("-i")){
                parameters_filename = "parameters.in";
                std::cout << "Using the input parameter file: " << parameters_filename << std::endl;
            }
            else {
                throw(0);
            }
        }
        else if (argc == 1){
            parameters_filename = "parameters.in";
            std::cout << "Using the input parameter file: " << parameters_filename << std::endl;
        }
        else {
            throw(1);
        }
    }
    catch (int e)
    {
        std::string e_mess;
        if (e == 0){e_mess = "Invalid command line option given. The only argument should be to specify the input file name.";}
        else if (e == 1){e_mess = "Too many command line arguments were given. The only argument should be to specify the input file name.";}
        else {e_mess = "Uknown exception.";}

        std::cerr << std::endl << std::endl
        << "----------------------------------------------------"
        << std::endl;
        std::cerr << "PRISMS-PF: Exception on processing: " << std::endl
        << e_mess << std::endl
        << "Aborting!" << std::endl
        << "----------------------------------------------------"
        << std::endl;
        return 1;
    }


    try
    {
        dealii::deallog.depth_console(0);

        // Before fully parsing the parameter file, we need to know how many field variables there are and whether they
        // are scalars or vectors, how many postprocessing variables there are, how many sets of elastic constants there are,
        // and how many user-defined constants there are.

        variableAttributeLoader variable_attributes;
        inputFileReader input_file_reader(parameters_filename,variable_attributes);

        // Continue based on the number of dimensions and degree of the elements specified in the input file
        switch (input_file_reader.number_of_dimensions)
        {
            case 2:
            {
                userInputParameters<2> userInputs(input_file_reader,input_file_reader.parameter_handler,variable_attributes);
                switch (userInputs.degree)
                {
                    case(1):
                    {
                        customPDE<2,1> problem(userInputs);
                        problem.buildFields();
                        problem.init ();
                        problem.solve();
                        break;
                    }
                    case(2):
                    {
                        customPDE<2,2> problem(userInputs);
                        problem.buildFields();
                        problem.init ();
                        problem.solve();
                        break;
                    }
                    case(3):
                    {
                        customPDE<2,3> problem(userInputs);
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
                userInputParameters<3> userInputs(input_file_reader,input_file_reader.parameter_handler,variable_attributes);
                switch (userInputs.degree)
                {
                    case(1):
                    {
                        customPDE<3,1> problem(userInputs);
                        problem.buildFields();
                        problem.init ();
                        problem.solve();
                        break;
                    }
                    case(2):
                    {
                        customPDE<3,2> problem(userInputs);
                        problem.buildFields();
                        problem.init ();
                        problem.solve();
                        break;
                    }
                    case(3):
                    {
                        customPDE<3,3> problem(userInputs);
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
        std::cerr << "PRISMS-PF: Exception on processing: " << std::endl
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
