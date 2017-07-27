// Methods for the userInputParameters class
#include "../../include/userInputParameters.h"
//#include "../../include/sortIndexEntryPairList.h"

template <int dim>
userInputParameters<dim>::userInputParameters(inputFileReader & input_file_reader, dealii::ParameterHandler & parameter_handler, variableAttributeLoader variable_attributes){
    loadVariableAttributes(variable_attributes);

    unsigned int _number_of_variables = input_file_reader.var_types.size();
    unsigned int _number_of_pp_variables = input_file_reader.num_pp_vars;

    // Load the inputs into the class member variables

    // Meshing parameters
    domain_size.push_back(parameter_handler.get_double("Domain size X"));
    if (dim > 1){
	       domain_size.push_back(parameter_handler.get_double("Domain size Y"));
           if (dim > 2){
               domain_size.push_back(parameter_handler.get_double("Domain size Z"));
           }
    }

    subdivisions.push_back(parameter_handler.get_integer("Subdivisions X"));
	if (dim > 1){
		subdivisions.push_back(parameter_handler.get_integer("Subdivisions Y"));
		if (dim > 2){
			subdivisions.push_back(parameter_handler.get_integer("Subdivisions Z"));
		}
	}

    refine_factor = parameter_handler.get_integer("Refine factor");

    degree = parameter_handler.get_integer("Element degree");

    // Adaptive meshing parameters
    h_adaptivity = parameter_handler.get_bool("Mesh adaptivity");
    max_refinement_level = parameter_handler.get_integer("Max refinement level");
    min_refinement_level = parameter_handler.get_integer("Min refinement level");

    // Use built-in deal.II utilities to split up a string and convert it to a vector of doubles or ints
    refine_criterion_fields = dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(parameter_handler.get("Refinement criteria fields")));
    refine_window_max = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Refinement window max")));
    refine_window_min = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Refinement window min")));

    skip_remeshing_steps = parameter_handler.get_integer("Steps between remeshing operations");

    // Time stepping parameters
    dtValue = parameter_handler.get_double("Time step");
    int totalIncrements_temp = parameter_handler.get_integer("Number of time steps");
    finalTime = parameter_handler.get_double("Simulation end time");

    // Elliptic solver parameters
    solver_type = parameter_handler.get("Linear solver");
    abs_tol = parameter_handler.get_bool("Use absolute convergence tolerance");
    solver_tolerance = parameter_handler.get_double("Solver tolerance value");
    max_solver_iterations = parameter_handler.get_integer("Maximum allowed solver iterations");

    // Output parameters
    std::string output_condition = parameter_handler.get("Output condition");
    unsigned int num_outputs = parameter_handler.get_integer("Number of outputs");
    std::vector<int> user_given_time_step_list_temp = dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(parameter_handler.get("List of time steps to output")));
    std::vector<unsigned int> user_given_time_step_list;
    for (unsigned int i=0; i<user_given_time_step_list_temp.size(); i++) user_given_time_step_list.push_back(user_given_time_step_list_temp[i]);

    skip_print_steps = parameter_handler.get_integer("Skip print steps");
    output_file_type = parameter_handler.get("Output file type");
    output_file_name = parameter_handler.get("Output file name (base)");
    calc_energy = parameter_handler.get_bool("Calculate the free energy");

    // Field variable definitions


    number_of_variables = _number_of_variables;


    // If all of the variables are ELLIPTIC, then totalIncrements should be 1 and finalTime should be 0
    bool only_elliptic_pdes = true;
    for (unsigned int i=0; i<var_eq_type.size(); i++){
        if (var_eq_type.at(i) == PARABOLIC){
            only_elliptic_pdes = false;
            break;
        }
    }

    // Determine the maximum number of time steps
    if (only_elliptic_pdes){
        totalIncrements = 1;
        finalTime = 0.0;
    }
    else{
        if ((totalIncrements_temp >= 0) && (finalTime >= 0.0)) {
            if (std::ceil(finalTime/dtValue) > totalIncrements_temp) {
                totalIncrements = totalIncrements_temp;
                finalTime = totalIncrements*dtValue;
            }
            else {
                totalIncrements = std::ceil(finalTime/dtValue);
            }
        }
        else if ((totalIncrements_temp >= 0) && (finalTime < 0.0)) {
            totalIncrements = totalIncrements_temp;
            finalTime = totalIncrements*dtValue;
        }
        else if ((totalIncrements_temp < 0) && (finalTime >= 0.0)) {
            totalIncrements = std::ceil(finalTime/dtValue);
        }
        else {
            // Should change to an exception
            std::cerr << "Invalid selections for the final time and the number of increments. At least one should be given in the input file and should be positive." << std::endl;
            std::cout << finalTime << " " << totalIncrements_temp << std::endl;
            abort();
        }
    }


    // Use these inputs to create a list of time steps where the code should output, stored in the member
    outputTimeStepList = setOutputTimeSteps(output_condition, num_outputs,user_given_time_step_list);

    // Variables for loading in PField ICs
    std::vector<std::string> load_ICs_temp = dealii::Utilities::split_string_list(parameter_handler.get("Load initial conditions"));
    std::vector<std::string> load_parallel_file_temp = dealii::Utilities::split_string_list(parameter_handler.get("Load parallel file"));

    if (boost::iequals(load_ICs_temp.at(0),"void")){
        for (unsigned int var=0; var<number_of_variables; var++){
            load_ICs.push_back(false);
            load_parallel_file.push_back(false);
        }
    }
    else {
        for (unsigned int var=0; var<number_of_variables; var++){
            if (boost::iequals(load_ICs_temp.at(var),"true")){
                load_ICs.push_back(true);
            }
            else {
                load_ICs.push_back(false);
            }
            if (boost::iequals(load_parallel_file_temp.at(var),"true")){
                load_parallel_file.push_back(true);
            }
            else {
                load_parallel_file.push_back(false);
            }
        }
    }

    std::vector<std::string> load_file_name = dealii::Utilities::split_string_list(parameter_handler.get("File names"));
    std::vector<std::string> load_field_name = dealii::Utilities::split_string_list(parameter_handler.get("Variable names in the files"));

    // Parameters for nucleation

    nucleus_semiaxes = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Nucleus semiaxes (x, y ,z)")));
    order_parameter_freeze_semiaxes = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Freeze zone semiaxes (x, y ,z)")));
    nucleus_hold_time = parameter_handler.get_double("Freeze time following nucleation");
    no_nucleation_border_thickness = parameter_handler.get_double("Nucleation-free border thickness");

    if (parameter_handler.get("Minimum allowed distance between nuclei") != ""){
        min_distance_between_nuclei = parameter_handler.get_double("Minimum allowed distance between nuclei");
    }
    else {
        min_distance_between_nuclei = 2.0 * (*(max_element(nucleus_semiaxes.begin(),nucleus_semiaxes.end())));
    }
    nucleation_order_parameter_cutoff = parameter_handler.get_double("Order parameter cutoff value");
    steps_between_nucleation_attempts = parameter_handler.get_integer("Time steps between nucleation attempts");


    // Load the boundary condition variables into list of BCs (where each element of the vector is one component of one variable)
    std::vector<std::string> list_of_BCs;
    for (unsigned int i=0; i<number_of_variables; i++){
        if (var_type[i] == SCALAR){
            std::string bc_text = "Boundary condition for variable ";
            bc_text.append(var_name.at(i));
            list_of_BCs.push_back(parameter_handler.get(bc_text));
        }
        else {
            std::string bc_text = "Boundary condition for variable ";
            bc_text.append(var_name.at(i));
            bc_text.append(", x component");
            list_of_BCs.push_back(parameter_handler.get(bc_text));

            bc_text = "Boundary condition for variable ";
            bc_text.append(var_name.at(i));
            bc_text.append(", y component");
            list_of_BCs.push_back(parameter_handler.get(bc_text));

            bc_text = "Boundary condition for variable ";
            bc_text.append(var_name.at(i));
            bc_text.append(", z component");
            list_of_BCs.push_back(parameter_handler.get(bc_text));
        }
    }

    // Load the BC information from the strings into a varBCs object
    load_BC_list(list_of_BCs);

    // Load the user-defined constants
    load_user_constants(input_file_reader,parameter_handler);

}


// Template instantiations
#include "../../include/userInputParameters_template_instantiations.h"
