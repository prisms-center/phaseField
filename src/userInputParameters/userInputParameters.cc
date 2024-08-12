// Methods for the userInputParameters class
#include "../../include/userInputParameters.h"
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

template <int dim>
userInputParameters<dim>::userInputParameters(inputFileReader & input_file_reader, dealii::ParameterHandler & parameter_handler, variableAttributeLoader variable_attributes){

    loadVariableAttributes(variable_attributes);

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
    skip_remeshing_steps = parameter_handler.get_integer("Steps between remeshing operations");

    max_refinement_level = parameter_handler.get_integer("Max refinement level");
    min_refinement_level = parameter_handler.get_integer("Min refinement level");

    // Enforce that the initial refinement level must be between the max and min level
    if (h_adaptivity && ((refine_factor < min_refinement_level) || (refine_factor > max_refinement_level))){
        std::cerr << "PRISMS-PF Error: The initial refinement factor must be between the maximum and minimum refinement levels when adaptive meshing is enabled." << std::endl;
        std::cerr << "Initial refinement level: " << refine_factor << " Maximum and minimum refinement levels: " << max_refinement_level << ", " << min_refinement_level << std::endl;
        abort();
    }

    // The adaptivity criterion for each variable has its own subsection
    for (unsigned int i=0; i<number_of_variables; i++){

        std::string subsection_text = "Refinement criterion: ";
        subsection_text.append(input_file_reader.var_names.at(i));

        parameter_handler.enter_subsection(subsection_text);
        {
            std::string crit_type_string = parameter_handler.get("Criterion type");
            if (crit_type_string.size() > 0){
                RefinementCriterion new_criterion;
                new_criterion.variable_index = i;
                new_criterion.variable_name = input_file_reader.var_names.at(i);
                if (boost::iequals(crit_type_string,"VALUE")){
                    new_criterion.criterion_type = VALUE;
                    new_criterion.value_lower_bound = parameter_handler.get_double("Value lower bound");
                    new_criterion.value_upper_bound = parameter_handler.get_double("Value upper bound");

                    // Check to make sure that the upper bound is greater than or equal to the lower bound
                    if (new_criterion.value_upper_bound < new_criterion.value_lower_bound){
                        std::cerr << "PRISMS-PF Error: The upper bound for refinement for variable " << new_criterion.variable_name << " is less than the lower bound. Please correct this in the parameters file." << std::endl;
                    }
                }
                else if (boost::iequals(crit_type_string,"GRADIENT")){
                    new_criterion.criterion_type = GRADIENT;
                    new_criterion.gradient_lower_bound = parameter_handler.get_double("Gradient magnitude lower bound");
                }
                else if (boost::iequals(crit_type_string,"VALUE_AND_GRADIENT")){
                    new_criterion.criterion_type = VALUE_AND_GRADIENT;
                    new_criterion.value_lower_bound = parameter_handler.get_double("Value lower bound");
                    new_criterion.value_upper_bound = parameter_handler.get_double("Value upper bound");
                    new_criterion.gradient_lower_bound = parameter_handler.get_double("Gradient magnitude lower bound");

                    // Check to make sure that the upper bound is greater than or equal to the lower bound
                    if (new_criterion.value_upper_bound < new_criterion.value_lower_bound){
                        std::cerr << "PRISMS-PF Error: The upper bound for refinement for variable " << new_criterion.variable_name << " is less than the lower bound. Please correct this in the parameters file." << std::endl;
                    }
                }
                else {
                    std::cerr << "PRISMS-PF Error: The refinement criteria type found in the parameters file, " << crit_type_string << ", is not an allowed type. The allowed types are VALUE, GRADIENT, VALUE_AND_GRADIENT" << std::endl;
                    abort();
                }
                refinement_criteria.push_back(new_criterion);
            }
        }
        parameter_handler.leave_subsection();
    }

    // Time stepping parameters
    dtValue = parameter_handler.get_double("Time step");
    int totalIncrements_temp = parameter_handler.get_integer("Number of time steps");
    finalTime = parameter_handler.get_double("Simulation end time");

    // Linear solver parameters
    for (unsigned int i=0; i<number_of_variables; i++){
        if (input_file_reader.var_eq_types.at(i) == TIME_INDEPENDENT || input_file_reader.var_eq_types.at(i) == IMPLICIT_TIME_DEPENDENT){
            std::string subsection_text = "Linear solver parameters: ";
            subsection_text.append(input_file_reader.var_names.at(i));

            parameter_handler.enter_subsection(subsection_text);
            {
                // Set the tolerance type
                SolverToleranceType temp_type;
                std::string type_string = parameter_handler.get("Tolerance type");
                if (boost::iequals(type_string,"ABSOLUTE_RESIDUAL")){
                    temp_type = ABSOLUTE_RESIDUAL;
                }
                else if (boost::iequals(type_string,"RELATIVE_RESIDUAL_CHANGE")){
                    temp_type = RELATIVE_RESIDUAL_CHANGE;
                }
                else if (boost::iequals(type_string,"ABSOLUTE_SOLUTION_CHANGE")){
                    temp_type = ABSOLUTE_SOLUTION_CHANGE;
                    std::cerr << "PRISMS-PF Error: Linear solver tolerance type " << type_string << " is not currently implemented, please use either ABSOLUTE_RESIDUAL or RELATIVE_RESIDUAL_CHANGE" << std::endl;
                    abort();
                }
                else {
                    std::cerr << "PRISMS-PF Error: Linear solver tolerance type " << type_string << " is not one of the allowed values (ABSOLUTE_RESIDUAL, RELATIVE_RESIDUAL_CHANGE, ABSOLUTE_SOLUTION_CHANGE)" << std::endl;
                    abort();
                }

                // Set the tolerance value
                double temp_value = parameter_handler.get_double("Tolerance value");

                // Set the maximum number of iterations
                unsigned int temp_max_iterations = parameter_handler.get_integer("Maximum linear solver iterations");

                linear_solver_parameters.loadParameters(i,temp_type,temp_value,temp_max_iterations);
            }
            parameter_handler.leave_subsection();
        }
    }


    // Non-linear solver parameters
    std::vector<bool> var_nonlinear = variable_attributes.var_nonlinear;

    nonlinear_solver_parameters.setMaxIterations(parameter_handler.get_integer("Maximum nonlinear solver iterations"));

    for (unsigned int i=0; i<var_nonlinear.size(); i++){
        if (var_nonlinear.at(i)){
            std::string subsection_text = "Nonlinear solver parameters: ";
            subsection_text.append(input_file_reader.var_names.at(i));

            parameter_handler.enter_subsection(subsection_text);
            {
                // Set the tolerance type
                SolverToleranceType temp_type;
                std::string type_string = parameter_handler.get("Tolerance type");
                if (boost::iequals(type_string,"ABSOLUTE_RESIDUAL")){
                    temp_type = ABSOLUTE_RESIDUAL;
                }
                else if (boost::iequals(type_string,"RELATIVE_RESIDUAL_CHANGE")){
                    temp_type = RELATIVE_RESIDUAL_CHANGE;
                }
                else if (boost::iequals(type_string,"ABSOLUTE_SOLUTION_CHANGE")){
                    temp_type = ABSOLUTE_SOLUTION_CHANGE;
                }
                else {
                    std::cerr << "PRISMS-PF Error: Nonlinear solver tolerance type " << type_string << " is not one of the allowed values (ABSOLUTE_RESIDUAL, RELATIVE_RESIDUAL_CHANGE, ABSOLUTE_SOLUTION_CHANGE)" << std::endl;
                    abort();
                }

                // Set the tolerance value
                double temp_value = parameter_handler.get_double("Tolerance value");

                // Set the backtrace damping flag
                bool temp_backtrack_damping = parameter_handler.get_bool("Use backtracking line search damping");

                // Set the backtracking step size modifier
                double temp_step_modifier = parameter_handler.get_double("Backtracking step size modifier");

                // Set the constant that determines how much the residual must decrease to be accepted as sufficient
                double temp_residual_decrease_coeff = parameter_handler.get_double("Backtracking residual decrease coefficient");

                // Set the default damping coefficient (used if backtracking isn't used)
                double temp_damping_coefficient = parameter_handler.get_double("Constant damping value");

                // Set whether to use the solution of Laplace's equation instead of the IC in ICs_and_BCs.h as the initial guess for nonlinear, time independent equations
                bool temp_laplace_for_initial_guess;
                if (var_eq_type[i] == TIME_INDEPENDENT){
                    temp_laplace_for_initial_guess = parameter_handler.get_bool("Use Laplace's equation to determine the initial guess");
                }
                else {
                    temp_laplace_for_initial_guess = false;
                    if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0){
                        std::cout << "PRISMS-PF Warning: Laplace's equation is only used to generate the initial guess for time independent equations. The equation for variable " << var_name[i] << " is not a time independent equation. No initial guess is needed for this equation." << std::endl;
                    }
                }

                nonlinear_solver_parameters.loadParameters(i,temp_type,temp_value,temp_backtrack_damping,temp_step_modifier,temp_residual_decrease_coeff,temp_damping_coefficient,temp_laplace_for_initial_guess);
            }
            parameter_handler.leave_subsection();
        }
    }

    // Set the max number of nonlinear iterations
    if (var_nonlinear.size() == 0){
        nonlinear_solver_parameters.setMaxIterations(0);
    }

    // Output parameters
    std::string output_condition = parameter_handler.get("Output condition");
    unsigned int num_outputs = parameter_handler.get_integer("Number of outputs");
    std::vector<int> user_given_time_step_list_temp = dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(parameter_handler.get("List of time steps to output")));
    std::vector<unsigned int> user_given_time_step_list;
    for (unsigned int i=0; i<user_given_time_step_list_temp.size(); i++) user_given_time_step_list.push_back(user_given_time_step_list_temp[i]);

    skip_print_steps = parameter_handler.get_integer("Skip print steps");
    output_file_type = parameter_handler.get("Output file type");
    output_file_name = parameter_handler.get("Output file name (base)");

    output_vtu_per_process = parameter_handler.get_bool("Output separate files per process");
    if ((output_file_type == "vtk") && (!output_vtu_per_process)){
        output_vtu_per_process = true;
        if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0){
            std::cout << "PRISMS-PF Warning: 'Output file type' given as 'vtk' and 'Output separate files per process' given as 'false'. Shared output files are not supported for the vtk output format. Separate files per process will be created." << std::endl;
        }
    }

    print_timing_with_output = parameter_handler.get_bool("Print timing information with output");

    // Field variable definitions

    // If all of the variables are ELLIPTIC, then totalIncrements should be 1 and finalTime should be 0
    bool only_time_independent_pdes = true;
    for (unsigned int i=0; i<var_eq_type.size(); i++){
        if (var_eq_type.at(i) == EXPLICIT_TIME_DEPENDENT || var_eq_type.at(i) == IMPLICIT_TIME_DEPENDENT){
            only_time_independent_pdes = false;
            break;
        }
    }

    // Determine the maximum number of time steps
    if (only_time_independent_pdes){
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

            parameter_handler.enter_subsection(nucleation_text);
            {
                unsigned int var_index = i;
                std::vector<double> semiaxes = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Nucleus semiaxes (x, y, z)")));
                std::vector<double> ellipsoid_rotation = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Nucleus rotation in degrees (x, y, z)")));
                std::vector<double> freeze_semiaxes = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Freeze zone semiaxes (x, y, z)")));
                double hold_time = parameter_handler.get_double("Freeze time following nucleation");
                double no_nucleation_border_thickness = parameter_handler.get_double("Nucleation-free border thickness");

                nucleationParameters<dim> temp(var_index,semiaxes,freeze_semiaxes,ellipsoid_rotation,hold_time,no_nucleation_border_thickness);
                nucleation_parameters_list.push_back(temp);

                // Validate nucleation input
                if (semiaxes.size() < dim || semiaxes.size() > 3){
                    std::cerr << "PRISMS-PF Error: The number of nucleus semiaxes given in the 'parameters.in' file must be at least the number of dimensions and no more than 3." << std::endl;
                    abort();
                }
                if (freeze_semiaxes.size() < dim || freeze_semiaxes.size() > 3){
                    std::cerr << "PRISMS-PF Error: The number of nucleation freeze zone semiaxes given in the 'parameters.in' file must be at least the number of dimensions and no more than 3." << std::endl;
                    abort();
                }
                if (ellipsoid_rotation.size() != 3){
                    std::cerr << "PRISMS-PF Error: Exactly three nucleus rotation angles must be given in the 'parameters.in' file." << std::endl;
                    abort();
                }
            }
            parameter_handler.leave_subsection();

        }
    }
    for (unsigned int i=0; i<nucleation_parameters_list.size(); i++){
        nucleation_parameters_list_index[nucleation_parameters_list.at(i).var_index] = i;
    }


    if (parameter_handler.get("Minimum allowed distance between nuclei") != "-1"){
        min_distance_between_nuclei = parameter_handler.get_double("Minimum allowed distance between nuclei");
    }
    else if (nucleation_parameters_list.size() > 1) {
        min_distance_between_nuclei = 2.0 * (*(max_element(nucleation_parameters_list[0].semiaxes.begin(),nucleation_parameters_list[0].semiaxes.end())));
    }
    min_distance_between_OP = parameter_handler.get_double("Minimum allowed distance between nuclei OP");
    multiple_nuclei_per_order_parameter = parameter_handler.get_bool("Allow multiple nuclei per order parameter");
    nucleation_order_parameter_cutoff = parameter_handler.get_double("Order parameter cutoff value");
    steps_between_nucleation_attempts = parameter_handler.get_integer("Time steps between nucleation attempts");
    nucleation_start_time = parameter_handler.get_double("Nucleation start time");
    nucleation_end_time = parameter_handler.get_double("Nucleation end time");

    // Load the grain remapping parameters
    grain_remapping_activated = parameter_handler.get_bool("Activate grain reassignment");

    skip_grain_reassignment_steps = parameter_handler.get_integer("Time steps between grain reassignments");

    order_parameter_threshold = parameter_handler.get_double("Order parameter cutoff for grain identification");

    buffer_between_grains = parameter_handler.get_double("Buffer between grains before reassignment");
    if (buffer_between_grains < 0.0 && grain_remapping_activated == true){
        std::cerr << "PRISMS-PF Error: If grain reassignment is activated, a non-negative buffer distance must be given. See the 'Buffer between grains before reassignment' entry in parameters.in." << std::endl;
        abort();
    }

    std::vector<std::string> variables_for_remapping_str = dealii::Utilities::split_string_list(parameter_handler.get("Order parameter fields for grain reassignment"));
    for (unsigned int field=0; field<variables_for_remapping_str.size(); field++){
        bool field_found = false;
        for (unsigned int i=0; i<number_of_variables; i++ ){
            if (boost::iequals(variables_for_remapping_str[field], variable_attributes.var_name_list[i].second)){
                variables_for_remapping.push_back(variable_attributes.var_name_list[i].first);
                field_found = true;
                break;
            }
        }
        if (field_found == false && grain_remapping_activated == true){
            std::cerr << "PRISMS-PF Error: Entries in the list of order parameter fields used for grain reassignment must match the variable names in equations.h." << std::endl;
            std::cerr << variables_for_remapping_str[field] << std::endl;
            abort();
        }
    }

    load_grain_structure = parameter_handler.get_bool("Load grain structure");
    load_unstructured_grid = parameter_handler.get_bool("Load as unstructured grid");
    grain_structure_filename = parameter_handler.get("Grain structure filename");
    grain_structure_variable_name = parameter_handler.get("Grain structure variable name");
    num_grain_smoothing_cycles = parameter_handler.get_integer("Number of smoothing cycles after grain structure loading");
    min_radius_for_loading_grains = parameter_handler.get_double("Minimum radius for loaded grains");

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
