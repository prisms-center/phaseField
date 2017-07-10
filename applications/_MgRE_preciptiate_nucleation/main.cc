// Coupled Cahn-Hilliard/Allen-Cahn example application with nucleation

// Header files
#include "../../include/initialConditions.h"
#include "../../include/initialCondition_template_instantiations.h"
#include "../../include/matrixFreePDE.h"
#include "../../src/models/mechanics/computeStress.h"
#include "../../include/inputFileReader.h"
#include "customPDE.h"
#include "equations.h"
#include "ICs_and_BCs.h"

// Header file for postprocessing that may or may not exist
#ifdef POSTPROCESS_FILE_EXISTS
#include "postprocess.h"
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
    try
    {
        dealii::deallog.depth_console(0);

        // Before fully parsing the parameter file, we need to know how many field variables there are and whether they
        // are scalars or vectors, how many postprocessing variables there are, how many sets of elastic constants there are,
        // and how many user-defined constants there are.
        inputFileReader input_file_reader;
        const std::vector<std::string> var_types = input_file_reader.get_subsection_entry_list("parameters.in","Variable","Variable type","SCALAR");
        const unsigned int num_materials = input_file_reader.get_number_of_entries("parameters.in","subsection","Material");
        const unsigned int num_pp_vars = input_file_reader.get_number_of_entries("parameters.in","subsection","Postprocessing variable");
        const unsigned int num_constants = input_file_reader.get_number_of_entries("parameters.in","set","Model constant");

        std::cout << "Number of model constants: " << num_constants << std::endl;

        // Read in all of the parameters now
        dealii::ParameterHandler parameter_handler;
        input_file_reader.declare_parameters(parameter_handler,var_types,
                                             num_materials,num_pp_vars,num_constants);
        parameter_handler.read_input("parameters.in");

        // Continue based on the number of dimensions and degree of the elements specified in the input file
        switch (parameter_handler.get_integer("Number of dimensions"))
        {
            case 2:
            {
                userInputParameters<2> userInputs;
                userInputs.loadInputParameters(parameter_handler,var_types.size(),num_materials,num_pp_vars,num_constants);
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
                userInputParameters<3> userInputs;
                userInputs.loadInputParameters(parameter_handler,var_types.size(),num_materials,num_pp_vars,num_constants);
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
