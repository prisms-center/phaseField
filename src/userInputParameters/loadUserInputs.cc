// Methods for the userInputParameters class

template <int dim>
void userInputParameters<dim>::loadUserInput(){

	// Meshing parameters
	domain_size.push_back(spanX);
	domain_size.push_back(spanY);
	domain_size.push_back(spanZ);

	subdivisions.push_back(subdivisionsX);
	if (dim > 1){
		subdivisions.push_back(subdivisionsY);
		if (dim > 2){
			subdivisions.push_back(subdivisionsZ);
		}
	}

	refine_factor = refineFactor;

	// Mesh refinement parameters
	h_adaptivity = hAdaptivity;
	max_refinement_level = maxRefinementLevel;
	min_refinement_level = minRefinementLevel;
	{int temp[] = refineCriterionFields;
	vectorLoad(temp, sizeof(temp), refine_criterion_fields);}
	{double temp[] = refineWindowMax;
	vectorLoad(temp, sizeof(temp), refine_window_max);}
	{double temp[] = refineWindowMin;
	vectorLoad(temp, sizeof(temp), refine_window_min);}

	skip_remeshing_steps = skipRemeshingSteps;

	// Output parameters
	write_output = writeOutput;
	output_condition = outputCondition;
	num_outputs = numOutputs;
	{unsigned int temp[] = outputList;
	vectorLoad(temp,sizeof(temp),user_given_time_step_list);}
	skip_print_steps = skipPrintSteps;
	output_file_type = outputFileType;

	calc_energy = calcEnergy;

	// Nucleation flag
	nucleation_occurs = nucleationOccurs;

	// Elliptic solver inputs
	solver_type = solverType;
	if (solver_type != "SolverCG"){
		std::cout << "Currently the only allowed solver is 'SolverCG'" << std::endl;
	}
	abs_tol = absTol;
	solver_tolerance = solverTolerance;
	max_solver_iterations = maxSolverIterations;

	//initialize time step variables
	#ifdef timeStep
	dtValue=timeStep;
	#endif
	#ifdef timeFinal
	finalTime=timeFinal;
	#else
	finalTime = 0.0;
	#endif


	// Determine the maximum number of time steps
	if (std::ceil(finalTime/timeStep) < timeIncrements){
		totalIncrements = std::ceil(finalTime/timeStep);
	}
	else {
		totalIncrements = timeIncrements;
	}

	// Load in inputs from equations.h

	// Somewhat convoluted initialization so as not to rely on C++11 initializer lists (which not all compiler have yet)
	{std::string temp_string[] = variable_name;
	vectorLoad(temp_string,sizeof(temp_string),var_name);}

	{std::string temp_string[] = variable_type;
	vectorLoad(temp_string,sizeof(temp_string),var_type);}

	{std::string temp_string[] = variable_eq_type;
	vectorLoad(temp_string,sizeof(temp_string),var_eq_type);}


	{bool temp[] = need_val;
	vectorLoad(temp, sizeof(temp), need_value);}

	{bool temp[] = need_grad;
	vectorLoad(temp, sizeof(temp), need_gradient);}

	{bool temp[] = need_hess;
	vectorLoad(temp, sizeof(temp), need_hessian);}

	{bool temp[] = need_val_residual;
	vectorLoad(temp, sizeof(temp), value_residual);}

	{bool temp[] = need_grad_residual;
	vectorLoad(temp, sizeof(temp), gradient_residual);}

	number_of_variables = var_name.size();


	#ifdef need_val_LHS
	{bool temp[] = need_val_LHS;
	vectorLoad(temp, sizeof(temp), need_value_LHS);}
	#else
	for (unsigned int i=0; i<number_of_variables; i++)
		need_value_LHS.push_back(false);
	#endif
	#ifdef need_grad_LHS
	{bool temp[] = need_grad_LHS;
	vectorLoad(temp, sizeof(temp), need_gradient_LHS);}
	#else
	for (unsigned int i=0; i<number_of_variables; i++)
		need_gradient_LHS.push_back(false);
	#endif
	#ifdef need_hess_LHS
	{bool temp[] = need_hess_LHS;
	vectorLoad(temp, sizeof(temp), need_hessian_LHS);}
	#else
	for (unsigned int i=0; i<number_of_variables; i++)
		need_hessian_LHS.push_back(false);
	#endif
	#ifdef need_val_residual_LHS
	{bool temp[] = need_val_residual_LHS;
	vectorLoad(temp, sizeof(temp), value_residual_LHS);}
	#else
	for (unsigned int i=0; i<number_of_variables; i++)
		value_residual_LHS.push_back(false);
	#endif
	#ifdef need_grad_residual_LHS
	{bool temp[] = need_grad_residual_LHS;
	vectorLoad(temp, sizeof(temp), gradient_residual_LHS);}
	#else
	for (unsigned int i=0; i<number_of_variables; i++)
		gradient_residual_LHS.push_back(false);
	#endif


	// initialize CIJ vector
	#if defined(MaterialModels) && defined(MaterialConstants)

	// Somewhat convoluted initialization so as not to rely on C++11 initializer lists (which not all compiler have yet)
	std::vector<std::string> temp_mat_models;
	{std::string temp_string[] = MaterialModels;
	vectorLoad(temp_string,sizeof(temp_string),temp_mat_models);}

	std::vector<std::vector<double> > temp_mat_consts;
	double temp_array[][21] = MaterialConstants; // Largest allowable length is 21, uninitialized slots are set to zero
	for (unsigned int num_mat=0; num_mat < temp_mat_models.size(); num_mat++){
		std::vector<double> temp_vec;
		vectorLoad(temp_array[num_mat], sizeof(temp_array[num_mat]), temp_vec);
		temp_mat_consts.push_back(temp_vec);
	}

	elasticityModel mat_model;

	dealii::ConditionalOStream pcout(std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0);

	dealii::Tensor<2, 2*dim-1+dim/3, dealii::VectorizedArray<double> > CIJ_temp;
	for (unsigned int mater_num=0; mater_num < temp_mat_consts.size(); mater_num++){
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
	#endif

	// If the LHS variable attributes aren't defined

	// If nucleation isn't specifically turned on, set nucleation_occurs to false
	#ifndef nucleation_occurs
		#define nucleation_occurs false
	#endif

	// Load variable information for calculating the RHS
	varInfoListRHS.reserve(number_of_variables);
	unsigned int scalar_var_index = 0;
	unsigned int vector_var_index = 0;
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

	// Initializing variables for loading in initial conditions using PFields
	#ifdef loadICs
		load_ICs = loadICs;
	#else
		for (unsigned int i=0; i<number_of_variables; i++){
			load_ICs.push_back(false);
		}
	#endif

	#ifdef loadSerialFile
		load_serial_file = loadSerialFile;
	#else
		for (unsigned int i=0; i<number_of_variables; i++){
			load_serial_file.push_back(false);
		}
	#endif
	#ifdef loadFileName
		load_file_name = loadFileName;
	#else
		for (unsigned int i=0; i<number_of_variables; i++){
			load_file_name.push_back("void");
		}
	#endif
	#ifdef loadFieldName
		load_field_name = loadFieldName;
	#else
		for (unsigned int i=0; i<number_of_variables; i++){
			load_field_name.push_back("void");
		}
	#endif

	// Postprocessing inputs

	// Somewhat convoluted initialization so as not to rely on C++11 initializer lists (which not all compilers have yet)
	#ifdef pp_variable_name

	{std::string temp_string[] = pp_variable_name;
	vectorLoad(temp_string,sizeof(temp_string),pp_var_name);}

	{std::string temp_string[] = pp_variable_type;
	vectorLoad(temp_string,sizeof(temp_string),pp_var_type);}

	{bool temp[] = pp_need_val;
	vectorLoad(temp, sizeof(temp), pp_need_value);}

	{bool temp[] = pp_need_grad;
	vectorLoad(temp, sizeof(temp), pp_need_gradient);}

	{bool temp[] = pp_need_hess;
	vectorLoad(temp, sizeof(temp), pp_need_hessian);}

	{bool temp[] = pp_need_val_residual;
	vectorLoad(temp, sizeof(temp), pp_value_residual);}

	{bool temp[] = pp_need_grad_residual;
	vectorLoad(temp, sizeof(temp), pp_gradient_residual);}

	pp_number_of_variables = pp_var_name.size();


	// Load variable information for calculating the RHS
	pp_varInfoList.reserve(pp_number_of_variables);
	scalar_var_index = 0;
	vector_var_index = 0;
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
	#endif

}




