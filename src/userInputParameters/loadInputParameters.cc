// Methods for the userInputParameters class
#include "../../include/userInputParameters.h"

template <int dim>
userInputParameters<dim>::userInputParameters(inputFileReader & input_file_reader, dealii::ParameterHandler & parameter_handler){
    unsigned int _number_of_variables = input_file_reader.var_types.size();
    unsigned int _number_of_materials = input_file_reader.num_materials;
    unsigned int _number_of_pp_variables = input_file_reader.num_pp_vars;
    unsigned int _number_of_constants = input_file_reader.num_constants;

    for (unsigned int i=0; i<input_file_reader.model_constant_names.size(); i++){
        model_constant_name_map[input_file_reader.model_constant_names[i]] = i;
    }

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
    std::vector<bool> nucleating_variable;
    std::vector<bool> need_value_nucleation;

    std::vector<bool> need_value;
	std::vector<bool> need_gradient;
	std::vector<bool> need_hessian;
	std::vector<bool> value_residual;
	std::vector<bool> gradient_residual;

    std::vector<bool> need_value_LHS;
	std::vector<bool> need_gradient_LHS;
	std::vector<bool> need_hessian_LHS;
	std::vector<bool> value_residual_LHS;
	std::vector<bool> gradient_residual_LHS;

    std::vector<bool> need_value_pp;
	std::vector<bool> need_gradient_pp;
	std::vector<bool> need_hessian_pp;

    number_of_variables = _number_of_variables;

    for (unsigned int i=0; i<number_of_variables; i++){
        std::string equation_text = "Variable ";
        equation_text.append(dealii::Utilities::int_to_string(i));

        parameter_handler.enter_subsection(equation_text);
        {
            var_name.push_back(parameter_handler.get("Variable name"));

            std::string var_type_temp = parameter_handler.get("Variable type");
            if (var_type_temp == "SCALAR"){
                var_type.push_back(SCALAR);
            }
            else if (var_type_temp == "VECTOR"){
                var_type.push_back(VECTOR);
            }
            else {
                std::cerr << "PRISMS-PF Error: Variable type must be 'SCALAR' or 'VECTOR'." << std::endl;
                abort();
            }

            std::string var_eq_type_temp = parameter_handler.get("Equation type");
            if (var_eq_type_temp == "PARABOLIC"){
                var_eq_type.push_back(PARABOLIC);
            }
            else if (var_eq_type_temp == "ELLIPTIC"){
                var_eq_type.push_back(ELLIPTIC);
            }
            else {
                std::cerr << "PRISMS-PF Error: Variable equation type must be 'PARABOLIC' or 'ELLIPTIC'." << std::endl;
                abort();
            }

            need_value.push_back(parameter_handler.get_bool("Need variable value"));
            need_gradient.push_back(parameter_handler.get_bool("Need variable gradient"));
            need_hessian.push_back(parameter_handler.get_bool("Need variable hessian"));
            value_residual.push_back(parameter_handler.get_bool("Need value residual term"));
            gradient_residual.push_back(parameter_handler.get_bool("Need gradient residual term"));

            need_value_LHS.push_back(parameter_handler.get_bool("Need variable value (LHS)"));
            need_gradient_LHS.push_back(parameter_handler.get_bool("Need variable gradient (LHS)"));
            need_hessian_LHS.push_back(parameter_handler.get_bool("Need variable hessian (LHS)"));
            value_residual_LHS.push_back(parameter_handler.get_bool("Need value residual term (LHS)"));
            gradient_residual_LHS.push_back(parameter_handler.get_bool("Need gradient residual term (LHS)"));

            need_value_pp.push_back(parameter_handler.get_bool("Need variable value (post-processing)"));
            need_gradient_pp.push_back(parameter_handler.get_bool("Need variable gradient (post-processing)"));
            need_hessian_pp.push_back(parameter_handler.get_bool("Need variable hessian (post-processing)"));

            nucleating_variable.push_back(parameter_handler.get_bool("Nucleating variable"));
            need_value_nucleation.push_back(parameter_handler.get_bool("Need variable value (nucleation)"));
        }
        parameter_handler.leave_subsection();
    }

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
    else {
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

    // Elastic constants
    unsigned int number_of_materials = _number_of_materials;

    std::vector<std::vector<double> > temp_mat_consts;
    std::vector<std::string> temp_mat_models;

    for (unsigned int i=0; i<number_of_materials; i++){
        std::string material_text = "Material ";
        material_text.append(dealii::Utilities::int_to_string(i));

        parameter_handler.enter_subsection(material_text);
        {
            std::vector<double> temp_vec;
            temp_vec = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic constants")));
            temp_mat_consts.push_back(temp_vec);
            temp_mat_models.push_back(parameter_handler.get("Material symmetry"));
        }
        parameter_handler.leave_subsection();
    }

    elasticityModel mat_model;

	dealii::ConditionalOStream pcout(std::cout, dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0);

	dealii::Tensor<2, 2*dim-1+dim/3, dealii::VectorizedArray<double> > CIJ_temp;
	for (unsigned int mater_num=0; mater_num < number_of_materials; mater_num++){
		if (temp_mat_models[mater_num] == "ISOTROPIC"){
			mat_model = ISOTROPIC;
		}
		else if (temp_mat_models[mater_num] == "TRANSVERSE"){
			mat_model = TRANSVERSE;
		}
		else if (temp_mat_models[mater_num] == "ORTHOTROPIC"){
			mat_model = ORTHOTROPIC;
		}
		else if (temp_mat_models[mater_num] == "ANISOTROPIC"){
			mat_model = ANISOTROPIC;
		}
		else {
			// Should change to an exception
			std::cerr << "Elastic material model is invalid, please use ISOTROPIC, TRANSVERSE, ORTHOTROPIC, or ANISOTROPIC" << std::endl;
		}

		getCIJMatrix<dim>(mat_model, temp_mat_consts[mater_num], CIJ_temp, pcout);
		CIJ_list.push_back(CIJ_temp);

		material_moduli.CIJ_list = CIJ_list;
	}

    // Variables for loading in PField ICs
    std::vector<std::string> load_ICs_temp = dealii::Utilities::split_string_list(parameter_handler.get("Load initial conditions"));
    std::vector<std::string> load_serial_file_temp = dealii::Utilities::split_string_list(parameter_handler.get("Load serial file"));

    if (boost::iequals(load_ICs_temp.at(0),"void")){
        for (unsigned int var=0; var<number_of_variables; var++){
            load_ICs.push_back(false);
            load_serial_file.push_back(false);
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
            if (boost::iequals(load_serial_file_temp.at(var),"true")){
                load_serial_file.push_back(true);
            }
            else {
                load_serial_file.push_back(false);
            }
        }
    }

    std::vector<std::string> load_file_name = dealii::Utilities::split_string_list(parameter_handler.get("File names"));
    std::vector<std::string> load_field_name = dealii::Utilities::split_string_list(parameter_handler.get("Variable names in the files"));

    // Parameters for postprocessing
    pp_number_of_variables = _number_of_pp_variables;

    std::vector<bool> pp_value_residual;
	std::vector<bool> pp_gradient_residual;


    if (pp_number_of_variables > 0){
        postProcessingRequired = true;
    }
    else {
        postProcessingRequired = false;
    }

    for (unsigned int i=0; i<pp_number_of_variables; i++){
        std::string pp_var_text = "Postprocessing variable ";
        pp_var_text.append(dealii::Utilities::int_to_string(i));

        parameter_handler.enter_subsection(pp_var_text);
        {
            pp_var_name.push_back(parameter_handler.get("Variable name"));

            std::string pp_var_type_temp = parameter_handler.get("Variable type");
            if (pp_var_type_temp == "SCALAR"){
                pp_var_type.push_back(SCALAR);
            }
            else if (pp_var_type_temp == "VECTOR"){
                pp_var_type_temp.push_back(VECTOR);
            }
            else {
                std::cerr << "PRISMS-PF Error: Variable type must be 'SCALAR' or 'VECTOR'." << std::endl;
                abort();
            }

            pp_value_residual.push_back(parameter_handler.get_bool("Need value residual term"));
            pp_gradient_residual.push_back(parameter_handler.get_bool("Need gradient residual term"));
        }
        parameter_handler.leave_subsection();
    }

    // Load variable information for postprocessing
    // First, the info list for the base field variables
	pp_baseVarInfoList.reserve(number_of_variables);
	unsigned int scalar_var_index = 0;
	unsigned int vector_var_index = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		variable_info varInfo;

        varInfo.need_value = need_value_pp[i];
        varInfo.need_gradient = need_gradient_pp[i];
        varInfo.need_hessian = need_hessian_pp[i];

        varInfo.global_var_index = i;

		if (need_value_pp[i] or need_gradient_pp[i] or need_hessian_pp[i]){
            varInfo.var_needed = true;
		}
        else {
            varInfo.var_needed = false;
        }

        if (var_type[i] == SCALAR){
            varInfo.is_scalar = true;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = scalar_var_index;
                scalar_var_index++;
            }
        }
        else {
            varInfo.is_scalar = false;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = vector_var_index;
                vector_var_index++;
            }
        }

        pp_baseVarInfoList.push_back(varInfo);
	}

    // Second, the info list for the postprocessing field variables
	pp_varInfoList.reserve(pp_number_of_variables);
	scalar_var_index = 0;
	vector_var_index = 0;
	for (unsigned int i=0; i<pp_number_of_variables; i++){
		variable_info varInfo;
        varInfo.var_needed = true;

        varInfo.value_residual = pp_value_residual[i];
        varInfo.gradient_residual = pp_gradient_residual[i];

		varInfo.global_var_index = i;
		if (pp_var_type[i] == SCALAR){
			varInfo.is_scalar = true;
			varInfo.scalar_or_vector_index = scalar_var_index;
			scalar_var_index++;
		}
		else {
			varInfo.is_scalar = false;
			varInfo.scalar_or_vector_index = vector_var_index;
			vector_var_index++;
		}
		pp_varInfoList.push_back(varInfo);
    }

    // Parameters for nucleation
    for (unsigned int i=0; i<number_of_variables; i++){
        if (nucleating_variable.at(i)==true){
            nucleating_variable_indices.push_back(i);
        }
        if (need_value_nucleation.at(i) || nucleating_variable.at(i)){
            nucleation_need_value.push_back(i);
        }
    }

    if (nucleating_variable_indices.size() > 0){
        nucleation_occurs = true;
    }
    else {
        nucleation_occurs = false;
    }

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
            bc_text.append(dealii::Utilities::int_to_string(i));
            list_of_BCs.push_back(parameter_handler.get(bc_text));
        }
        else {
            std::string bc_text = "Boundary condition for variable ";
            bc_text.append(dealii::Utilities::int_to_string(i));
            bc_text.append(", x component");
            list_of_BCs.push_back(parameter_handler.get(bc_text));

            bc_text = "Boundary condition for variable ";
            bc_text.append(dealii::Utilities::int_to_string(i));
            bc_text.append(", y component");
            list_of_BCs.push_back(parameter_handler.get(bc_text));

            bc_text = "Boundary condition for variable ";
            bc_text.append(dealii::Utilities::int_to_string(i));
            bc_text.append(", z component");
            list_of_BCs.push_back(parameter_handler.get(bc_text));
        }
    }

    // Load the BC information from the strings into a varBCs object
    load_BC_list(list_of_BCs);

    // Load the user-defined constants
    for (unsigned int i=0; i<_number_of_constants; i++){
        std::string constants_text = "Model constant ";
        //constants_text.append(dealii::Utilities::int_to_string(i));
        constants_text.append(input_file_reader.model_constant_names[i]);
        std::vector<std::string> model_constants_strings = dealii::Utilities::split_string_list(parameter_handler.get(constants_text));
        if (model_constants_strings.size() == 2){
            if (boost::iequals(model_constants_strings.at(1),"double")){
                model_constants.push_back(dealii::Utilities::string_to_double(model_constants_strings.at(0)));

            }
            else if (boost::iequals(model_constants_strings.at(1),"int")){
                model_constants.push_back(dealii::Utilities::string_to_int(model_constants_strings.at(0)));
            }
            else if (boost::iequals(model_constants_strings.at(1),"bool")){
                bool temp;
                if (boost::iequals(model_constants_strings.at(0),"true")){
                    temp = true;
                }
                else {
                    temp = false;
                }
                model_constants.push_back(temp);
            }
            else {
                std::cerr << "PRISMS-PF ERROR: The type for user-defined variables must be 'double', 'int', 'bool', or 'tensor'." << std::endl;
                abort();
            }
        }
        else if (model_constants_strings.size() < 2){
            std::cerr << "PRISMS-PF ERROR: Users must input two fields for user-defined variables (value and type)." << std::endl;
            abort();
        }
        else {
            if (boost::iequals(model_constants_strings.at(model_constants_strings.size()-1),"tensor")){
                unsigned int num_elements = model_constants_strings.size()-1;

                // Strip parentheses from the input, counting how many rows there are
                unsigned int open_parentheses = 0;
                unsigned int close_parentheses = 0;
                for (unsigned int element=0; element<num_elements; element++){
                    std::size_t index = 0;
                    while (index != std::string::npos){
                        index = model_constants_strings.at(element).find("(");
                        if (index != std::string::npos) {
                            model_constants_strings.at(element).erase(index,1);
                            open_parentheses++;
                        }
                    }
                    index = 0;
                    while (index != std::string::npos){
                        index = model_constants_strings.at(element).find(")");
                        if (index != std::string::npos) {
                            model_constants_strings.at(element).erase(index,1);
                            close_parentheses++;
                        }
                    }
                }
                if (open_parentheses != close_parentheses){
                    std::cerr << "PRISMS-PF ERROR: User-defined constant tensor does not have the same number of open and close parentheses." << std::endl;
                    abort();
                }
                // Rank 1 tensor
                if (open_parentheses < 3){
                    if (num_elements > 1 && num_elements < 4){
                        dealii::Tensor<1,dim> temp;
                        for (unsigned int i=0; i<dim; i++){
                            temp[i] = dealii::Utilities::string_to_double(model_constants_strings.at(i));
                        }
                        model_constants.push_back(temp);
                    }
                    else {
                        std::cerr << "PRISMS-PF ERROR: The columns in user-defined constant tensors cannot be longer than 3 elements (internally truncated to the number of dimensions)." << std::endl;
                        abort();
                    }
                }
                // Rank 2 tensor
                else if (open_parentheses < 5){
                    unsigned int row_length;
                    if (num_elements == 4){
                        row_length = 2;
                        if (dim > 2){
                            std::cerr << "PRISMS-PF ERROR: User-defined constant tensor does not have enough elements, for 3D calculations matrices must be 3x3." << std::endl;
                            abort();
                        }
                    }
                    else if (num_elements == 9){
                        row_length = 3;
                    }
                    else {
                        std::cerr << "PRISMS-PF ERROR: User-defined constant tensor does not have the correct number of elements, matrices must be 2x2 or 3x3." << std::endl;
                        abort();
                    }

                    dealii::Tensor<2,dim> temp;
                    for (unsigned int i=0; i<dim; i++){
                        for (unsigned int j=0; j<dim; j++){
                            temp[i][j] = dealii::Utilities::string_to_double(model_constants_strings.at(i*row_length+j));
                        }
                    }
                    model_constants.push_back(temp);
                }
            }
            else {
                std::cerr << "PRISMS-PF ERROR: Only user-defined constant tensors may have multiple elements." << std::endl;
                abort();
            }
        }
    }

    // --------------------------------------------------------------------------------------------------------
    // Build the varInfoList objects using the parameters loaded from the input file
    // --------------------------------------------------------------------------------------------------------

    // Load variable information for calculating the RHS
	varInfoListRHS.reserve(number_of_variables);
	scalar_var_index = 0;
	vector_var_index = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		variable_info varInfo;

        varInfo.need_value = need_value[i];
        varInfo.need_gradient = need_gradient[i];
        varInfo.need_hessian = need_hessian[i];
        varInfo.value_residual = value_residual[i];
        varInfo.gradient_residual = gradient_residual[i];

        varInfo.global_var_index = i;

		if (need_value[i] or need_gradient[i] or need_hessian[i]){
            varInfo.var_needed = true;
		}
        else {
            varInfo.var_needed = false;
        }

        if (var_type[i] == SCALAR){
            varInfo.is_scalar = true;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = scalar_var_index;
                scalar_var_index++;
            }
        }
        else {
            varInfo.is_scalar = false;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = vector_var_index;
                vector_var_index++;
            }
        }

        varInfoListRHS.push_back(varInfo);
	}

	// Load variable information for calculating the LHS
	num_var_LHS = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		if (need_value_LHS[i] or need_gradient_LHS[i] or need_hessian_LHS[i]){
			num_var_LHS++;
		}
	}

	varInfoListLHS.reserve(num_var_LHS);
	scalar_var_index = 0;
	vector_var_index = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		variable_info varInfo;

        varInfo.need_value = need_value_LHS[i];
        varInfo.need_gradient = need_gradient_LHS[i];
        varInfo.need_hessian = need_hessian_LHS[i];
        varInfo.value_residual = value_residual_LHS[i];
        varInfo.gradient_residual = gradient_residual_LHS[i];

        varInfo.global_var_index = i;


		if (need_value_LHS[i] or need_gradient_LHS[i] or need_hessian_LHS[i]){
            varInfo.var_needed = true;
		}
        else {
            varInfo.var_needed = false;
		}

        if (var_type[i] == SCALAR){
            varInfo.is_scalar = true;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = scalar_var_index;
                scalar_var_index++;
            }
        }
        else {
            varInfo.is_scalar = false;
            if (varInfo.var_needed){
                varInfo.scalar_or_vector_index = vector_var_index;
                vector_var_index++;
            }
        }

        varInfoListLHS.push_back(varInfo);
	}

}

template <int dim>
void userInputParameters<dim>::load_BC_list(std::vector<std::string> list_of_BCs){
    // Load the BC information from the strings into a varBCs object
    // Move this to a new method and write a unit test for it!!!!

    std::vector<std::string> temp;

    for (unsigned int i=0; i<list_of_BCs.size(); i++){
        varBCs<dim> newBC;
        temp = dealii::Utilities::split_string_list(list_of_BCs[i]);

        // If there is only one BC listed, make another dim*2-1 copies of it so that the same BC is applied for all boundaries
        if (temp.size() == 1){
            for (unsigned int boundary=0; boundary<(dim*2-1); boundary++){
                temp.push_back(temp[0]);
            }
        }

        // Load the BC for each boundary into 'newBC'
        for (unsigned int i=0; i<(2*dim); i++){
            if (boost::iequals(temp[i],"ZERO_DERIVATIVE")){
                newBC.var_BC_type.push_back(ZERO_DERIVATIVE);
                newBC.var_BC_val.push_back(0.0);
            }
            else if (boost::iequals(temp[i],"PERIODIC")){
                newBC.var_BC_type.push_back(PERIODIC);
                newBC.var_BC_val.push_back(0.0);
            }
            else if (boost::iequals(temp[i].substr(0,9),"DIRICHLET")){
                newBC.var_BC_type.push_back(DIRICHLET);
                std::string dirichlet_val = temp[i].substr(10,temp[i].size());
                dirichlet_val = dealii::Utilities::trim(dirichlet_val);
                newBC.var_BC_val.push_back(dealii::Utilities::string_to_double(dirichlet_val));
            }
            else {
                std::cout << temp[i].substr(0,8) << std::endl;
                std::cout << "Error: Boundary conditions specified improperly." << std::endl;
                abort();
            }
        }
        BC_list.push_back(newBC);

            // Validate input using something like this:
            // try{
            //     if ((BC_type_dim1_min == "PERIODIC") && (BC_type_dim1_max != "PERIODIC")){
            //         throw 0;
            //     }
            //     if ((BC_type_dim2_min == "PERIODIC") && (BC_type_dim2_max != "PERIODIC")){
            //         throw 0;
            //     }
            //     if ((BC_type_dim3_min == "PERIODIC") && (BC_type_dim3_max != "PERIODIC")){
            //         throw 0;
            //     }
            // }
            // catch (int e){
            //     if (e == 0){
            //         std::cout << "Error: For periodic BCs, both faces for a given direction must be set as periodic. "
            //                 "Please check the BCs that are set in ICs_and_BCs.h." << std::endl;
            //     }
            //     abort();
            // }
    }

}

template <int dim>
std::vector<unsigned int> userInputParameters<dim>::setOutputTimeSteps(const std::string outputSpacingType, unsigned int numberOfOutputs,
                    const std::vector<unsigned int> & userGivenTimeStepList)
{
    std::vector<unsigned int> timeStepList;

    if (numberOfOutputs > 0) {
        if (outputSpacingType == "EQUAL_SPACING"){
            if (numberOfOutputs > totalIncrements)
            numberOfOutputs = totalIncrements;

            for (unsigned int iter = 0; iter <= totalIncrements; iter += totalIncrements/numberOfOutputs){
                timeStepList.push_back(iter);
            }
        }
        else if (outputSpacingType == "LOG_SPACING"){
            timeStepList.push_back(0);
            for (unsigned int output = 1; output <= numberOfOutputs; output++){
                timeStepList.push_back(round(std::pow(10,double(output)/double(numberOfOutputs)*std::log10(totalIncrements))));
            }
        }
        else if (outputSpacingType == "N_PER_DECADE"){
            timeStepList.push_back(0);
            timeStepList.push_back(1);
            for (unsigned int iter = 2; iter <= totalIncrements; iter++){
                int decade = std::ceil(std::log10(iter));
                int step_size = (std::pow(10,decade))/numberOfOutputs;
                if (iter%step_size == 0){
                    timeStepList.push_back(iter);
                }

            }
        }
        else if (outputSpacingType == "LIST"){
            timeStepList = userGivenTimeStepList;
        }
    }
    else {
        timeStepList.push_back(totalIncrements+1);
    }

    return timeStepList;
}


// Template instantiations
template class userInputParameters<2>;
template class userInputParameters<3>;
