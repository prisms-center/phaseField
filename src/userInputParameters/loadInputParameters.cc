// Methods for the userInputParameters class
#include "../../include/userInputParameters.h"

template <int dim>
void userInputParameters<dim>::loadInputParameters(dealii::ParameterHandler & parameter_handler,
                                                    unsigned int _number_of_variables, unsigned int _number_of_materials, unsigned int _number_of_pp_variables){

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

    // Nucleation parameters
    nucleation_occurs = parameter_handler.get_bool("Allow nucleation");

    // Field variable definitions
    number_of_variables = _number_of_variables;

    for (unsigned int i=0; i<number_of_variables; i++){
        std::string equation_text = "Equation ";
        equation_text.append(dealii::Utilities::int_to_string(i));

        parameter_handler.enter_subsection(equation_text);
        {
            var_name.push_back(parameter_handler.get("Variable name"));
            var_type.push_back(parameter_handler.get("Variable type"));
            var_eq_type.push_back(parameter_handler.get("Equation type"));

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
        }
        parameter_handler.leave_subsection();
    }

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
			std::cout << "Elastic material model is invalid, please use ISOTROPIC, TRANSVERSE, ORTHOTROPIC, or ANISOTROPIC" << std::endl;
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
            pp_var_type.push_back(parameter_handler.get("Variable type"));

            pp_need_value.push_back(parameter_handler.get_bool("Need variable value"));
            pp_need_gradient.push_back(parameter_handler.get_bool("Need variable gradient"));
            pp_need_hessian.push_back(parameter_handler.get_bool("Need variable hessian"));
            pp_value_residual.push_back(parameter_handler.get_bool("Need value residual term"));
            pp_gradient_residual.push_back(parameter_handler.get_bool("Need gradient residual term"));
        }
        parameter_handler.leave_subsection();
    }


    // Load variable information for calculating the RHS
	pp_varInfoList.reserve(pp_number_of_variables);
	unsigned int scalar_var_index = 0;
	unsigned int vector_var_index = 0;
	for (unsigned int i=0; i<pp_number_of_variables; i++){
		variable_info varInfo;

		varInfo.global_var_index = i;
		if (pp_var_type[i] == "SCALAR"){
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

    // Load the boundary condition variables into list of BCs (where each element of the vector is one component of one variable)
    std::vector<std::string> list_of_BCs;
    for (unsigned int i=0; i<number_of_variables; i++){
        if (boost::iequals(var_type[i],"SCALAR")){
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
    load_BC_list(list_of_BCs, BC_list);

    // --------------------------------------------------------------------------------------------------------
    // Build the varInfoList objects using the parameters loaded from the input file
    // --------------------------------------------------------------------------------------------------------

    // Load variable information for calculating the RHS
	varInfoListRHS.reserve(number_of_variables);
	scalar_var_index = 0;
	vector_var_index = 0;
	for (unsigned int i=0; i<number_of_variables; i++){
		variable_info varInfo;
		if (need_value[i] or need_gradient[i] or need_hessian[i]){
			varInfo.global_var_index = i;
			if (var_type[i] == "SCALAR"){
				varInfo.is_scalar = true;
				varInfo.scalar_or_vector_index = scalar_var_index;
				scalar_var_index++;
			}
			else {
				varInfo.is_scalar = false;
				varInfo.scalar_or_vector_index = vector_var_index;
				vector_var_index++;
			}
			varInfoListRHS.push_back(varInfo);
		}
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
		if (need_value_LHS[i] or need_gradient_LHS[i] or need_hessian_LHS[i]){
			varInfo.global_var_index = i;
			if (var_type[i] == "SCALAR"){
				varInfo.is_scalar = true;
				varInfo.scalar_or_vector_index = scalar_var_index;
				scalar_var_index++;
			}
			else {
				varInfo.is_scalar = false;
				varInfo.scalar_or_vector_index = vector_var_index;
				vector_var_index++;
			}
			varInfoListLHS.push_back(varInfo);
		}
	}

}

template <int dim>
void userInputParameters<dim>::load_BC_list(std::vector<std::string> list_of_BCs, std::vector<varBCs<dim> > & BC_list){
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

// Template instantiations
template class userInputParameters<2>;
template class userInputParameters<3>;
