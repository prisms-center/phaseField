// Allen-Cahn example application

// Header files
#include "../../include/initialConditions.h"
#include "../../include/matrixFreePDE.h"
#include "customPDE.h"
#include "equations.h"
#include "ICs_and_BCs.h"
#include "postprocess.h"
#include "../../include/initialCondition_template_instantiations.h"
#include "../../include/inputFileReader.h"

//main
int main (int argc, char **argv)
{
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,dealii::numbers::invalid_unsigned_int);
    try
    {
        dealii::deallog.depth_console(0);

        // Before fully parsing the parameter file, we need to know how many field variables there are and whether they
        // are scalars or vectors, how many postprocessing variables there are, and how many sets of elastic constants there are.
        inputFileReader input_file_reader;
        const std::vector<std::string> var_types = input_file_reader.get_subsection_entry_list("parameters.in","Equation","Variable type","SCALAR");
        const std::vector<std::string> material_types = input_file_reader.get_subsection_entry_list("parameters.in","Material","Material symmetry","ISOTROPIC");
        const std::vector<std::string> pp_var_types = input_file_reader.get_subsection_entry_list("parameters.in","Postprocessing variable","Variable type","SCALAR");

        // Read in all of the parameters now
        dealii::ParameterHandler parameter_handler;
        input_file_reader.declare_parameters(parameter_handler,"parameters.in",var_types,
                                             material_types.size(),pp_var_types.size());
        parameter_handler.read_input("parameters.in");

        // Continue based on the number of dimensions and degree of the elements specified in the input file
        switch (parameter_handler.get_integer("Number of dimensions"))
        {
            case 2:
            {
                userInputParameters<2> userInputs;
                userInputs.loadInputParameters(parameter_handler,var_types.size(),material_types.size(),pp_var_types.size());
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
                userInputs.loadInputParameters(parameter_handler,var_types.size(),material_types.size(),pp_var_types.size());
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
