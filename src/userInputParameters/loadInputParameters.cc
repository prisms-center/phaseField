// Methods for the userInputParameters class

template <int dim>
void userInputParameters<dim>::loadInputParameters(dealii::ParameterHandler & parameter_handler){

    // Load the inputs into the class member variables

    // Meshing parameters
    domain_size.push_back(parameter_handler.get_integer("Domain size X"));
    if (dim > 1){
	       domain_size.push_back(parameter_handler.get_integer("Domain size Y"));
           if (dim > 2){
               domain_size.push_back(parameter_handler.get_integer("Domain size Z"));
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

    // Determine the maximum number of time steps
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
        std::cout << "Invalid selections for the final time and the number of increments. At least one should be given in the input file and should be positive." << std::endl;
        abort();
    }

    // Elliptic solver parameters
    solver_type = parameter_handler.get("Linear solver");
    abs_tol = parameter_handler.get_bool("Use absolute convergence tolerance");
    solver_tolerance = parameter_handler.get_double("Solver tolerance value");
    max_solver_iterations = parameter_handler.get_integer("Maximum allowed solver iterations");

    // Output parameters
    write_output = parameter_handler.get_bool("Write output");
    output_condition = parameter_handler.get("Output condition");
    num_outputs = parameter_handler.get_integer("Number of outputs");

    std::vector<int> user_given_time_step_list_temp = dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(parameter_handler.get("List of time steps to output")));
    for (unsigned int i=0; i<user_given_time_step_list_temp.size(); i++) user_given_time_step_list.push_back(user_given_time_step_list_temp[i]);

    skip_print_steps = parameter_handler.get_integer("Skip print steps");
    output_file_type = parameter_handler.get("Output file type");
    calc_energy = parameter_handler.get_bool("Calculate the free energy");

    // Nucleation
    nucleation_occurs = parameter_handler.get_bool("Allow nucleation");

}
