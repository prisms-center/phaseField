// Methods for the userInputParameters class
#include "../../include/userInputParameters.h"
#include "../../include/sortIndexEntryPairList.h"

template <int dim>
userInputParameters<dim>::userInputParameters(inputFileReader & input_file_reader, dealii::ParameterHandler & parameter_handler, variableAttributeLoader variable_attributes){
    loadVariableAttributes(variable_attributes);

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
    std::vector<std::string> load_serial_file_temp = dealii::Utilities::split_string_list(parameter_handler.get("Load parallel file"));

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
        constants_text.append(input_file_reader.model_constant_names[i]);
        //std::cout << input_file_reader.model_constant_names[i] << std::endl;
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

template <int dim>
void userInputParameters<dim>::loadVariableAttributes(variableAttributeLoader variable_attributes){
    number_of_variables = variable_attributes.var_name_list.size();

    var_name = sortIndexEntryPairList(variable_attributes.var_name_list,number_of_variables,"var");
    var_type = sortIndexEntryPairList(variable_attributes.var_type_list,number_of_variables,SCALAR);
    var_eq_type = sortIndexEntryPairList(variable_attributes.var_eq_type_list,number_of_variables,PARABOLIC);

    // Sort the variable attributes and load them into individual vectors
    std::vector<bool> nucleating_variable = sortIndexEntryPairList(variable_attributes.nucleating_variable_list,number_of_variables,false);
    std::vector<bool> need_value_nucleation = sortIndexEntryPairList(variable_attributes.need_value_list_nucleation,number_of_variables,false);

    std::vector<bool> need_value = sortIndexEntryPairList(variable_attributes.need_value_list,number_of_variables,true);
	std::vector<bool> need_gradient = sortIndexEntryPairList(variable_attributes.need_gradient_list,number_of_variables,true);
	std::vector<bool> need_hessian = sortIndexEntryPairList(variable_attributes.need_hessian_list,number_of_variables,false);
	std::vector<bool> value_residual = sortIndexEntryPairList(variable_attributes.need_value_residual_list,number_of_variables,false);
	std::vector<bool> gradient_residual = sortIndexEntryPairList(variable_attributes.need_gradient_residual_list,number_of_variables,false);

    std::vector<bool> need_value_LHS = sortIndexEntryPairList(variable_attributes.need_value_list_LHS,number_of_variables,false);
	std::vector<bool> need_gradient_LHS = sortIndexEntryPairList(variable_attributes.need_gradient_list_LHS,number_of_variables,false);
	std::vector<bool> need_hessian_LHS = sortIndexEntryPairList(variable_attributes.need_hessian_list_LHS,number_of_variables,false);
	std::vector<bool> value_residual_LHS = sortIndexEntryPairList(variable_attributes.need_value_residual_list_LHS,number_of_variables,false);
	std::vector<bool> gradient_residual_LHS = sortIndexEntryPairList(variable_attributes.need_gradient_residual_list_LHS,number_of_variables,false);

    std::vector<bool> need_value_pp = sortIndexEntryPairList(variable_attributes.need_value_list_PP,number_of_variables,true);
	std::vector<bool> need_gradient_pp = sortIndexEntryPairList(variable_attributes.need_gradient_list_PP,number_of_variables,true);
	std::vector<bool> need_hessian_pp = sortIndexEntryPairList(variable_attributes.need_hessian_list_PP,number_of_variables,true);

    // Sort the post-processing variable attributes
    pp_number_of_variables = variable_attributes.var_name_list_PP.size();
    pp_var_name = sortIndexEntryPairList(variable_attributes.var_name_list_PP,pp_number_of_variables,"var");
    pp_var_type = sortIndexEntryPairList(variable_attributes.var_type_list_PP,pp_number_of_variables,SCALAR);
    std::vector<bool> pp_value_residual = sortIndexEntryPairList(variable_attributes.need_value_residual_list_PP,pp_number_of_variables,false);
	std::vector<bool> pp_gradient_residual = sortIndexEntryPairList(variable_attributes.need_gradient_residual_list_PP,pp_number_of_variables,false);
    pp_calc_integral = sortIndexEntryPairList(variable_attributes.output_integral_list,pp_number_of_variables,false);


    // Load some nucleation parameters
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


    // Load these attributes into the varInfoList objects
	varInfoListRHS.reserve(number_of_variables);
	unsigned int scalar_var_index = 0;
	unsigned int vector_var_index = 0;
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

    // Load variable information for postprocessing
    // First, the info list for the base field variables
	pp_baseVarInfoList.reserve(number_of_variables);
	scalar_var_index = 0;
	vector_var_index = 0;
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

    // Now load the information for the post-processing variables
    // Parameters for postprocessing
    if (pp_number_of_variables > 0){
        postProcessingRequired = true;
    }
    else {
        postProcessingRequired = false;
    }

    num_integrated_fields = 0;
    for (unsigned int i=0; i<pp_number_of_variables; i++){
        if (pp_calc_integral[i]){
            num_integrated_fields++;
            integrated_field_indices.push_back(i);
        }
    }

    // The info list for the postprocessing field variables
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

}


// Template instantiations
template class userInputParameters<2>;
template class userInputParameters<3>;
