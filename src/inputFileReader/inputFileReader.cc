// Methods for the inputFileReader class
#include "../../include/inputFileReader.h"

// Constructor
inputFileReader::inputFileReader(std::string input_file_name){
    var_types = get_subsection_entry_list("parameters.in","Variable","Variable type","SCALAR");
    num_materials = get_number_of_entries("parameters.in","subsection","Material");
    num_pp_vars = get_number_of_entries("parameters.in","subsection","Postprocessing variable");
    num_constants = get_number_of_entries("parameters.in","set","Model constant");

    model_constant_names = get_entry_name_ending_list("parameters.in","set", "Model constant");

    std::cout << "Number of constants: " << num_constants << std::endl;

    // Read in all of the parameters now
    declare_parameters(parameter_handler,var_types,
                                         num_materials,num_pp_vars,num_constants);
    parameter_handler.read_input("parameters.in");
    number_of_dimensions = parameter_handler.get_integer("Number of dimensions");
}


// Method to parse a single line to find a target key value pair
bool inputFileReader::parse_line(std::string line, const std::string keyword, const std::string entry_name,
                                    std::string & out_string, const bool expect_equals_sign) const
    {

    // Strip spaces at the front and back
    while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
    line.erase(0, 1);
    while ((line.size() > 0)
    && (line[line.size() - 1] == ' ' || line[line.size() - 1] == '\t'))
    line.erase(line.size() - 1, std::string::npos);

    // now see whether the line starts with 'keyword' followed by multiple spaces
    // if not, try next line (if the entry is "", then zero spaces after the keyword is ok)
    if (line.size() < keyword.size())
    return false;

    for (unsigned int i=0; i<keyword.size(); i++){
        if (line[i] != keyword[i])
            return false;
    }
    if (entry_name.size() > 0){
        if (!(line[keyword.size()] == ' ' || line[keyword.size()] == '\t'))
            return false;
    }

    // delete the "keyword" and then delete more spaces if present
    line.erase(0, keyword.size());
    while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
    line.erase(0, 1);
    // now see whether the next word is the word we look for
    if (line.find(entry_name) != 0)
    return false;

    line.erase(0, entry_name.size());
    while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
    line.erase(0, 1);

    // we'd expect an equals size here if expect_equals_sign is true
    if (expect_equals_sign){
        if ((line.size() < 1) || (line[0] != '='))
        return false;
    }

    // remove comment
    std::string::size_type pos = line.find('#');
    if (pos != std::string::npos)
    line.erase (pos);

    // trim the equals sign at the beginning and possibly following spaces
    // as well as spaces at the end
    if (expect_equals_sign)
        line.erase(0, 1);

    while ((line.size() > 0) && (line[0] == ' ' || line[0] == '\t'))
    line.erase(0, 1);

    while ((line.size() > 0) && (line[line.size()-1] == ' ' || line[line.size()-1] == '\t'))
    line.erase(line.size()-1, std::string::npos);

    out_string = line;
    return true;
}

// Method to parse an input file to get a list of variables from related subsections
std::vector<std::string> inputFileReader::get_subsection_entry_list(const std::string parameters_file_name,
                                                                    const std::string subsec_name, const std::string entry_name, const std::string default_entry) const
    {

    std::ifstream input_file;
    input_file.open(parameters_file_name);

    std::string line, entry;
    bool in_subsection = false;
    bool found_entry, desired_entry_found;
    unsigned int subsection_index;
    std::vector<std::string> entry_list;
    std::vector<unsigned int> index_list;

    // Loop through each line
    while (std::getline(input_file, line))
    {
        // If the line is the start of a subsection, turn 'in_subsection' to true and store the subsection index
        if (!in_subsection){
            found_entry = parse_line(line, "subsection",subsec_name, entry, false);
            if ( found_entry){
                in_subsection = true;
                subsection_index = dealii::Utilities::string_to_int(entry);
                desired_entry_found = false;
            }
        }
        // If in a subsection, look for the line setting the entry or for the end of the subsection
        else {
            found_entry = parse_line(line, "set", entry_name, entry, true);
            if (found_entry) {
                entry_list.push_back(entry);
                index_list.push_back(subsection_index);
                desired_entry_found = true;
            }
            found_entry = parse_line(line, "end", "", entry, false);
            if (found_entry) {
                if (!desired_entry_found){
                    entry_list.push_back(default_entry);
                    index_list.push_back(subsection_index);
                }
                in_subsection = false;
                desired_entry_found = false;
            }
        }
    }

    // Now sort the entry list vector so that it is in index order
    std::vector<std::string> sorted_entry_list;
    for (unsigned int i=0; i<entry_list.size(); i++){
        for (unsigned int j=0; j<entry_list.size(); j++){
            if (i == j){
                sorted_entry_list.push_back(entry_list[index_list[j]]);
                break;
            }
        }
    }
    return sorted_entry_list;
}

// Method to parse an input file to get a list of variables from related subsections
unsigned int inputFileReader::get_number_of_entries(const std::string parameters_file_name,
                                                                    const std::string keyword, const std::string entry_name) const
    {

    std::ifstream input_file;
    input_file.open(parameters_file_name);

    std::string line, entry;
    bool found_entry;

    unsigned int count = 0;

    // Loop through each line
    while (std::getline(input_file, line))
    {

        found_entry = parse_line(line, keyword, entry_name, entry, false);
        if (found_entry)
            count++;
    }
    return count;
}

// Method to parse an input file to get a list of variables from related subsections
std::vector<std::string> inputFileReader::get_entry_name_ending_list(const std::string parameters_file_name,
                                                                    const std::string keyword, const std::string entry_name_begining) const
    {

    std::ifstream input_file;
    input_file.open(parameters_file_name);

    std::string line, entry;
    bool found_entry;

    std::vector<std::string> entry_name_end_list;

    // Loop through each line
    while (std::getline(input_file, line))
    {

        found_entry = parse_line(line, keyword, entry_name_begining, entry, false);
        if (found_entry){

            // Strip whitespace, the equals sign, and everything after the equals sign

            // Strip whitespace at the beginning
            while ((entry.size() > 0) && (entry[0] == ' ' || entry[0] == '\t'))
            entry.erase(0, 1);

            // Strip everything up to the equals sign
            while ((entry.size() > 0) && (entry[entry.size() - 1] != '=' ))
            entry.erase(entry.size() - 1, std::string::npos);

            //Strip the equals sign
            entry.erase(entry.size() - 1, std::string::npos);

            // Strip whitespace between the entry name and the equals sign
            while ((entry.size() > 0) && (entry[entry.size() - 1] == ' ' || entry[entry.size() - 1] == '\t'))
            entry.erase(entry.size() - 1, std::string::npos);

            // Add it to the list
            entry_name_end_list.push_back(entry);
        }

    }
    return entry_name_end_list;
}

// Before fully parsing the parameter file, we need to know how many field variables there are
// This function is largely taken from ASPECT (https://github.com/geodynamics/aspect/blob/master/source/main.cc)
std::string inputFileReader::get_last_value_of_parameter(const std::string &parameters, const std::string &parameter_name) const
    {
        std::string return_value;

        //std::istringstream x_file(parameters);

        std::ifstream input_file;
        input_file.open(parameters);

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



void inputFileReader::declare_parameters(dealii::ParameterHandler & parameter_handler,
                                            const std::vector<std::string> var_types, const unsigned int number_of_materials,
                                            const unsigned int number_of_pp_variables, const unsigned int num_of_constants) const {

    // Declare all of the entries
    parameter_handler.declare_entry("Number of dimensions","-1",dealii::Patterns::Integer(),"The number of dimensions for the simulation.");

    parameter_handler.declare_entry("Domain size X","-1",dealii::Patterns::Double(),"The size of the domain in the x direction.");
    parameter_handler.declare_entry("Domain size Y","-1",dealii::Patterns::Double(),"The size of the domain in the y direction.");
    parameter_handler.declare_entry("Domain size Z","-1",dealii::Patterns::Double(),"The size of the domain in the z direction.");
    parameter_handler.declare_entry("Element degree","1",dealii::Patterns::Integer(),"The polynomial order of the finte element.");
    parameter_handler.declare_entry("Subdivisions X","1",dealii::Patterns::Integer(),"The number of mesh subdivisions in the x direction.");
    parameter_handler.declare_entry("Subdivisions Y","1",dealii::Patterns::Integer(),"The number of mesh subdivisions in the y direction.");
    parameter_handler.declare_entry("Subdivisions Z","1",dealii::Patterns::Integer(),"The number of mesh subdivisions in the z direction.");
    parameter_handler.declare_entry("Refine factor","-1",dealii::Patterns::Integer(),"The number of initial refinements of the coarse mesh.");

    parameter_handler.declare_entry("Mesh adaptivity","false",dealii::Patterns::Bool(),"Whether to enable mesh adaptivity.");
    parameter_handler.declare_entry("Max refinement level","-1",dealii::Patterns::Integer(),"The maximum level of refinement.");
    parameter_handler.declare_entry("Min refinement level","-1",dealii::Patterns::Integer(),"The minimum level of refinement.");
    parameter_handler.declare_entry("Refinement criteria fields","0",dealii::Patterns::List(dealii::Patterns::Integer()),"The list of fields used to determine mesh refinement.");
    parameter_handler.declare_entry("Refinement window max","",dealii::Patterns::List(dealii::Patterns::Anything()),"The upper limit for refinement for each of the criteria fields.");
    parameter_handler.declare_entry("Refinement window min","",dealii::Patterns::List(dealii::Patterns::Anything()),"The lower limit for refinement for each of the criteria fields.");
    parameter_handler.declare_entry("Steps between remeshing operations","1",dealii::Patterns::Integer(),"The number of time steps between mesh refinement operations.");

    parameter_handler.declare_entry("Number of time steps","-1",dealii::Patterns::Integer(),"The time step size for the simulation.");
    parameter_handler.declare_entry("Time step","-0.1",dealii::Patterns::Double(),"The time step size for the simulation.");
    parameter_handler.declare_entry("Simulation end time","-0.1",dealii::Patterns::Double(),"The value of simulated time where the simulation ends.");

    parameter_handler.declare_entry("Linear solver","SolverCG",dealii::Patterns::Anything(),"The linear solver (currently only SolverCG).");
    parameter_handler.declare_entry("Use absolute convergence tolerance","false",dealii::Patterns::Bool(),"Whether to use an absolute tolerance for the linear solver (versus a relative tolerance).");
    parameter_handler.declare_entry("Solver tolerance value","1.0e-3",dealii::Patterns::Double(),"The tolerance for the linear solver (either absolute or relative).");
    parameter_handler.declare_entry("Maximum allowed solver iterations","10000",dealii::Patterns::Integer(),"The maximum allowed number of iterations the linear solver is given to converge before being forced to exit.");

    parameter_handler.declare_entry("Output file name (base)","solution",dealii::Patterns::Anything(),"The name for the output file, before the time step and processor info are added.");
    parameter_handler.declare_entry("Output file type","vtu",dealii::Patterns::Anything(),"The output file type (either vtu or vtk).");
    parameter_handler.declare_entry("Output condition","EQUAL_SPACING",dealii::Patterns::Anything(),"The time step size for the simulation.");
    parameter_handler.declare_entry("List of time steps to output","0",dealii::Patterns::Anything(),"The list of time steps to output, used for the LIST type.");
    parameter_handler.declare_entry("Number of outputs","10",dealii::Patterns::Integer(),"The number of outputs (or number of outputs per decade for the N_PER_DECADE type).");
    parameter_handler.declare_entry("Skip print steps","1",dealii::Patterns::Integer(),"The number of time steps between updates to the screen.");
    parameter_handler.declare_entry("Calculate the free energy","false",dealii::Patterns::Bool(),"Whether to calculate the integrated free energy.");

    parameter_handler.declare_entry("Allow nucleation","false",dealii::Patterns::Bool(),"Whether to enable the explicit nucleation capabilties.");

    // Declare entries regarding the governing equations
    for (unsigned int i=0; i<var_types.size(); i++){
        std::string var_text = "Variable ";
        var_text.append(dealii::Utilities::int_to_string(i));

        parameter_handler.enter_subsection(var_text);
        {
            parameter_handler.declare_entry("Variable name","var",dealii::Patterns::Anything(),"The name of the field variable.");
            parameter_handler.declare_entry("Variable type","SCALAR",dealii::Patterns::Anything(),"Whether the variable is a SCALAR or a VECTOR.");
            parameter_handler.declare_entry("Equation type","PARABOLIC",dealii::Patterns::Anything(),"Whether the governing equation is PARABOLIC or ELLIPTIC.");

            parameter_handler.declare_entry("Need variable value","true",dealii::Patterns::Bool(),"Whether the value of the variable is needed for any of the residual equations.");
            parameter_handler.declare_entry("Need variable gradient","true",dealii::Patterns::Bool(),"Whether the gradient of the variable is needed for any of the residual equations.");
            parameter_handler.declare_entry("Need variable hessian","false",dealii::Patterns::Bool(),"Whether the hessian of the variable is needed for any of the residual equations.");

            parameter_handler.declare_entry("Need value residual term","true",dealii::Patterns::Bool(),"Whether the residual equation has a term proportional to the value of the test function.");
            parameter_handler.declare_entry("Need gradient residual term","true",dealii::Patterns::Bool(),"Whether the residual equation has a term proportional to the gradient of the test function.");

            parameter_handler.declare_entry("Need variable value (LHS)","false",dealii::Patterns::Bool(),"Whether the value of the variable is needed for any of the LHS residual equations.");
            parameter_handler.declare_entry("Need variable gradient (LHS)","false",dealii::Patterns::Bool(),"Whether the gradient of the variable is needed for any of the LHS residual equations.");
            parameter_handler.declare_entry("Need variable hessian (LHS)","false",dealii::Patterns::Bool(),"Whether the hessian of the variable is needed for any of the LHS residual equations.");

            parameter_handler.declare_entry("Need value residual term (LHS)","false",dealii::Patterns::Bool(),"Whether the LHS residual equation has a term proportional to the value of the test function.");
            parameter_handler.declare_entry("Need gradient residual term (LHS)","false",dealii::Patterns::Bool(),"Whether the LHS residual equation has a term proportional to the gradient of the test function.");

            parameter_handler.declare_entry("Need variable value (post-processing)","true",dealii::Patterns::Bool(),"Whether the value of the variable is needed for any of the post-processing residual equations.");
            parameter_handler.declare_entry("Need variable gradient (post-processing)","true",dealii::Patterns::Bool(),"Whether the gradient of the variable is needed for any of the post-processing residual equations.");
            parameter_handler.declare_entry("Need variable hessian (post-processing)","true",dealii::Patterns::Bool(),"Whether the hessian of the variable is needed for any of the post-processing residual equations.");

            parameter_handler.declare_entry("Nucleating variable","false",dealii::Patterns::Bool(),"Whether the variable is an order parameter that is allowed to nucleation.");
            parameter_handler.declare_entry("Need variable value (nucleation)","false",dealii::Patterns::Bool(),"Whether the value of the variable is needed for nucleation calculations (assumed to be true if the variable is a nucleating order parameter).");

        }
        parameter_handler.leave_subsection();
    }

    // Declare entries regarding the elastic constants
    for (unsigned int i=0; i<number_of_materials; i++){
        std::string material_text = "Material ";
        material_text.append(dealii::Utilities::int_to_string(i));

        parameter_handler.enter_subsection(material_text);
        {
            parameter_handler.declare_entry("Elastic constants","-1.0",dealii::Patterns::List(dealii::Patterns::Double()),"The elastic constants for the material (see the user guide for the ordering).");
            parameter_handler.declare_entry("Material symmetry","ISOTROPIC",dealii::Patterns::Anything(),"The elastic symmetry of the material (see the user guide for the list of options).");
        }
        parameter_handler.leave_subsection();
    }

    // Declare entries for reading initial conditions from file
    parameter_handler.declare_entry("Load initial conditions","void",dealii::Patterns::Anything(),"Whether to load the initial conditions for each variable from file.");
    parameter_handler.declare_entry("Load parallel file","void",dealii::Patterns::Anything(),"Whether all processors should read from a single file (versus each reading from separate files).");
    parameter_handler.declare_entry("File names","void",dealii::Patterns::Anything(),"The file name to load from for each variable.");
    parameter_handler.declare_entry("Variable names in the files","void",dealii::Patterns::Anything(),"What each variable is named in the field being loaded.");

    // Declare entries for the postprocessing variables
    for (unsigned int i=0; i<number_of_pp_variables; i++){
        std::string pp_var_text = "Postprocessing variable ";
        pp_var_text.append(dealii::Utilities::int_to_string(i));

        parameter_handler.enter_subsection(pp_var_text);
        {
            parameter_handler.declare_entry("Variable name","var",dealii::Patterns::Anything(),"The name of the field variable.");
            parameter_handler.declare_entry("Variable type","SCALAR",dealii::Patterns::Anything(),"Whether the variable is a SCALAR or a VECTOR.");

            parameter_handler.declare_entry("Need value residual term","true",dealii::Patterns::Bool(),"Whether the residual equation has a term proportional to the value of the test function.");
            parameter_handler.declare_entry("Need gradient residual term","true",dealii::Patterns::Bool(),"Whether the residual equation has a term proportional to the gradient of the test function.");

            parameter_handler.declare_entry("Output integral of the variable","false",dealii::Patterns::Bool(),"Whether to output the integral of the post-processed variable.");
        }
        parameter_handler.leave_subsection();
    }

    // Declare the boundary condition variables
    for (unsigned int i=0; i<var_types.size(); i++){
        if (boost::iequals(var_types[i],"SCALAR")){
            std::string bc_text = "Boundary condition for variable ";
            bc_text.append(dealii::Utilities::int_to_string(i));
            parameter_handler.declare_entry(bc_text,"",dealii::Patterns::Anything(),"The boundary conditions for one of the governing equations).");
        }
        else {
            std::string bc_text = "Boundary condition for variable ";
            bc_text.append(dealii::Utilities::int_to_string(i));
            bc_text.append(", x component");
            parameter_handler.declare_entry(bc_text,"",dealii::Patterns::Anything(),"The boundary conditions for one of the governing equations).");

            bc_text = "Boundary condition for variable ";
            bc_text.append(dealii::Utilities::int_to_string(i));
            bc_text.append(", y component");
            parameter_handler.declare_entry(bc_text,"",dealii::Patterns::Anything(),"The boundary conditions for one of the governing equations).");

            bc_text = "Boundary condition for variable ";
            bc_text.append(dealii::Utilities::int_to_string(i));
            bc_text.append(", z component");
            parameter_handler.declare_entry(bc_text,"",dealii::Patterns::Anything(),"The boundary conditions for one of the governing equations).");
        }

    }

    // Declare the nucleation parameters
    parameter_handler.declare_entry("Nucleus semiaxes (x, y ,z)","",dealii::Patterns::List(dealii::Patterns::Double()),"The semiaxes for nuclei placed with the explicit nucleation algorithm).");
    parameter_handler.declare_entry("Freeze zone semiaxes (x, y ,z)","",dealii::Patterns::List(dealii::Patterns::Double()),"The semiaxes for region where the order parameter is frozen for a period of time after placement.");
    parameter_handler.declare_entry("Freeze time following nucleation","0.0",dealii::Patterns::Double(),"Duration that the order parameter is frozen after placement.");
    parameter_handler.declare_entry("Nucleation-free border thickness","0.0",dealii::Patterns::Double(),"The thickness of the nucleation-free region near the domain boundaries (ignored for periodic BCs).");
    parameter_handler.declare_entry("Minimum allowed distance between nuclei","-1.0",dealii::Patterns::Double(),"The minimum allowed distance between nuclei placed during the same time step.");
    parameter_handler.declare_entry("Order parameter cutoff value","0.01",dealii::Patterns::Double(),"Order parameter cutoff value for nucleation (when the sum of all order parameters is above this value, no nucleation is attempted).");
    parameter_handler.declare_entry("Time steps between nucleation attempts","100",dealii::Patterns::Integer(),"The number of time steps between nucleation attempts.");

    // Declare the user-defined constants
    for (unsigned int i=0; i<num_of_constants; i++){
        std::string constants_text = "Model constant ";
        //constants_text.append(dealii::Utilities::int_to_string(i));
        constants_text.append(model_constant_names[i]);
        parameter_handler.declare_entry(constants_text,"0",dealii::Patterns::Anything(),"The value of a user-defined constant.");
    }



}
