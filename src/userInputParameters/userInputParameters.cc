// Methods for the userInputParameters class
#include "../../include/userInputParameters.h"
//#include "../../include/sortIndexEntryPairList.h"

template <int dim>
userInputParameters<dim>::userInputParameters(inputFileReader & input_file_reader, dealii::ParameterHandler & parameter_handler, variableAttributeLoader variable_attributes){
    loadVariableAttributes(variable_attributes);

    unsigned int _number_of_variables = input_file_reader.var_types.size();

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

    // Enforce that the initial refinement level must be between the max and min level
    if (h_adaptivity && ((refine_factor < min_refinement_level) || (refine_factor > max_refinement_level))){
        std::cerr << "PRISMS-PF Error: The initial refinement factor must be between the maximum and minimum refinement levels when adaptive meshing is enabled." << std::endl;
        std::cerr << "Initial refinement level: " << refine_factor << " Maximum and minimum refinement levels: " << max_refinement_level << ", " << min_refinement_level << std::endl;
        abort();
    }

    // Use built-in deal.II utilities to split up a string and convert it to a vector of doubles or ints
    std::vector<std::string> refine_criterion_fields_str = dealii::Utilities::split_string_list(parameter_handler.get("Refinement criteria fields"));
    for (unsigned int ref_field=0; ref_field<refine_criterion_fields_str.size(); ref_field++){
        bool field_found = false;
        for (unsigned int i=0; i<_number_of_variables; i++ ){
            if (boost::iequals(refine_criterion_fields_str[ref_field], variable_attributes.var_name_list[i].second)){
                refine_criterion_fields.push_back(variable_attributes.var_name_list[i].first);
                field_found = true;
                break;
            }
        }
        if (field_found == false && h_adaptivity == true){
            std::cerr << "PRISMS-PF Error: Entries in the list of fields used for refinement must match the variable names in equations.h." << std::endl;
            std::cerr << refine_criterion_fields_str[ref_field] << std::endl;
            abort();
        }
    }
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
    outputTimeStepList = setTimeStepList(output_condition, num_outputs,user_given_time_step_list);

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

    load_file_name = dealii::Utilities::split_string_list(parameter_handler.get("File names"));
    load_field_name = dealii::Utilities::split_string_list(parameter_handler.get("Variable names in the files"));

    // Parameters for checkpoint/restart
    resume_from_checkpoint = parameter_handler.get_bool("Load from a checkpoint");
    std::string checkpoint_condition = parameter_handler.get("Checkpoint condition");
    unsigned int num_checkpoints = parameter_handler.get_integer("Number of checkpoints");

    std::vector<int> user_given_checkpoint_time_step_list_temp = dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(parameter_handler.get("List of time steps to save checkpoints")));
    std::vector<unsigned int> user_given_checkpoint_time_step_list;
    for (unsigned int i=0; i<user_given_checkpoint_time_step_list_temp.size(); i++) user_given_checkpoint_time_step_list.push_back(user_given_checkpoint_time_step_list_temp[i]);

    checkpointTimeStepList = setTimeStepList(checkpoint_condition, num_checkpoints,user_given_checkpoint_time_step_list);

    // Parameters for nucleation

    for (unsigned int i=0; i<input_file_reader.var_types.size(); i++){
        if (input_file_reader.var_nucleates.at(i)){
            std::string nucleation_text = "Nucleation parameters: ";
            nucleation_text.append(input_file_reader.var_names.at(i));
            nucleationParameters<dim> temp;
            parameter_handler.enter_subsection(nucleation_text);
            {
                temp.var_index = i;
                temp.semiaxes = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Nucleus semiaxes (x, y, z)")));
                temp.ellipsoid_rotation = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Nucleus rotation in degrees (x, y, z)")));
                temp.freeze_semiaxes = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Freeze zone semiaxes (x, y, z)")));
                temp.hold_time = parameter_handler.get_double("Freeze time following nucleation");
                temp.no_nucleation_border_thickness = parameter_handler.get_double("Nucleation-free border thickness");

                temp.set_rotation_matrix();
            }
            parameter_handler.leave_subsection();
            nucleation_parameters_list.push_back(temp);
        }
    }

    if (parameter_handler.get("Minimum allowed distance between nuclei") != "-1"){
        min_distance_between_nuclei = parameter_handler.get_double("Minimum allowed distance between nuclei");
    }
    else if (nucleation_parameters_list.size() > 1) {
        min_distance_between_nuclei = 2.0 * (*(max_element(nucleation_parameters_list[0].semiaxes.begin(),nucleation_parameters_list[0].semiaxes.end())));
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
