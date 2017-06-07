// Allen-Cahn example application

// Header files
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

// Before fully parsing the parameter file, we need to know how many field variables there are
// This function is largely taken from ASPECT (https://github.com/geodynamics/aspect/blob/master/source/main.cc)
std::string
get_last_value_of_parameter(const std::string &parameters,
    const std::string &parameter_name)
    {
        std::string return_value;

        //std::istringstream x_file(parameters);

        std::ifstream input_file;
        input_file.open("parameters.in");
        
        std::string line;
        while (std::getline(input_file, line))
        {
            // get one line and strip spaces at the front and back
            while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
            line.erase(0, 1);
            while ((line.size() > 0)
            && (line[line.size() - 1] == ' ' || line[line.size() - 1] == '\t'))
            line.erase(line.size() - 1, std::string::npos);
            // now see whether the line starts with 'set' followed by multiple spaces
            // if not, try next line
            if (line.size() < 4)
            continue;

            if ((line[0] != 's') || (line[1] != 'e') || (line[2] != 't')
            || !(line[3] == ' ' || line[3] == '\t'))
            continue;

            // delete the "set " and then delete more spaces if present
            line.erase(0, 4);
            while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
            line.erase(0, 1);
            // now see whether the next word is the word we look for
            if (line.find(parameter_name) != 0)
            continue;

            line.erase(0, parameter_name.size());
            while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
            line.erase(0, 1);

            // we'd expect an equals size here
            if ((line.size() < 1) || (line[0] != '='))
            continue;

            // remove comment
            std::string::size_type pos = line.find('#');
            if (pos != std::string::npos)
            line.erase (pos);

            // trim the equals sign at the beginning and possibly following spaces
            // as well as spaces at the end
            line.erase(0, 1);
            while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
            line.erase(0, 1);
            while ((line.size() > 0) && (line[line.size()-1] == ' ' || line[line.size()-1] == '\t'))
            line.erase(line.size()-1, std::string::npos);

            // the rest should now be what we were looking for
            return_value = line;
        }

        input_file.close();

        return return_value;
    }



//main
int main (int argc, char **argv)
{
    dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,dealii::numbers::invalid_unsigned_int);
    try
    {
        dealii::deallog.depth_console(0);

        // Before fully parsing the parameter file, we need to know how many field variables there are, how many
        // postprocessing variables there are, and how many sets of elastic constants there are.
        const std::string numVar = get_last_value_of_parameter("parameters.in", "Number of equations");
        const std::string numMat = get_last_value_of_parameter("parameters.in", "Number of materials");
        const std::string numVar_pp = get_last_value_of_parameter("parameters.in", "Number of postprocessing variables");

        // Read in the parameters
        dealii::ParameterHandler parameter_handler;
        inputFileReader input_file_reader(parameter_handler,"parameters.in",dealii::Utilities::string_to_int(numVar),
                                            dealii::Utilities::string_to_int(numMat),dealii::Utilities::string_to_int(numVar_pp));

        const unsigned int dimension = parameter_handler.get_integer("Number of dimensions");

        switch (dimension)
        {
            case 2:
            {
                userInputParameters<2> userInputs;
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
