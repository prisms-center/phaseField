// =====================================================================
// CONSTRUCTOR
// =====================================================================

template <int dim>
generalizedProblem<dim>::generalizedProblem(): MatrixFreePDE<dim>()
{

#ifndef	timeIncrements
#define timeIncrements 1
#endif
#ifndef	timeStep
#define timeStep 0.0
#endif

	  //initialize time step variables
	#ifdef timeStep
	this->dtValue=timeStep;
	#endif
	#ifdef timeFinal
	this->finalTime=timeFinal;
	#else
	this->finalTime = 0.0;
	#endif


	  // Determine the maximum number of time steps
	if (std::ceil(this->finalTime/timeStep) < timeIncrements){
		this->totalIncrements = std::ceil(this->finalTime/timeStep);
	}
	else {
		this->totalIncrements = timeIncrements;
		}


	var_name = variable_name;
	var_type = variable_type;
	var_eq_type = variable_eq_type;

	need_value = need_val;
	need_gradient = need_grad;
	need_hessian = need_hess;
	value_residual = need_val_residual;
	gradient_residual = need_grad_residual;

	#ifdef need_val_LHS
	need_value_LHS = need_val_LHS;
	#else
	for (unsigned int i=0; i<num_var; i++)
		need_value_LHS.push_back(false);
	#endif
	#ifdef need_grad_LHS
	need_gradient_LHS = need_grad_LHS;
	#else
	for (unsigned int i=0; i<num_var; i++)
		need_gradient_LHS.push_back(false);
	#endif
	#ifdef need_hess_LHS
	need_hessian_LHS = need_hess_LHS;
	#else
	for (unsigned int i=0; i<num_var; i++)
		need_hessian_LHS.push_back(false);
	#endif
	#ifdef need_val_residual_LHS
	value_residual_LHS = need_val_residual_LHS;
	#else
	for (unsigned int i=0; i<num_var; i++)
		value_residual_LHS.push_back(false);
	#endif
	#ifdef need_grad_residual_LHS
	gradient_residual_LHS = need_grad_residual_LHS;
	#else
	for (unsigned int i=0; i<num_var; i++)
		gradient_residual_LHS.push_back(false);
	#endif


// initialize CIJ vector
#if defined(MaterialModels) && defined(MaterialConstants)
	std::vector<std::vector<double> > temp_mat_consts = MaterialConstants;
	std::vector<std::string> temp_mat_models = MaterialModels;
	elasticityModel mat_model;

	dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ_temp;
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

		getCIJMatrix<dim>(mat_model, temp_mat_consts[mater_num], CIJ_temp, this->pcout);
		CIJ_list.push_back(CIJ_temp);
	}
#endif


// I should probably get rid of this or move it, since it is only relevant to the precipitate case
c_dependent_misfit = false;
#if defined(sfts_linear1) && defined(sfts_linear2) && defined(sfts_linear3)
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		if ((std::abs(sfts_linear1[i][j])>1.0e-12)||(std::abs(sfts_linear2[i][j])>1.0e-12)||(std::abs(sfts_linear3[i][j])>1.0e-12)){
			c_dependent_misfit = true;
		}
	}
}
#endif


// If the LHS variable attributes aren't defined

// If nucleation isn't specifically turned on, set nucleation_occurs to false
#ifndef nucleation_occurs
	#define nucleation_occurs false
#endif

// Load variable information for calculating the RHS
varInfoListRHS.reserve(num_var);
unsigned int scalar_var_index = 0;
unsigned int vector_var_index = 0;
for (unsigned int i=0; i<num_var; i++){
	variable_info<dim> varInfo;
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
for (unsigned int i=0; i<num_var; i++){
	if (need_value_LHS[i] or need_gradient_LHS[i] or need_hessian_LHS[i]){
		num_var_LHS++;
	}
}

varInfoListLHS.reserve(num_var_LHS);
scalar_var_index = 0;
vector_var_index = 0;
for (unsigned int i=0; i<num_var; i++){
	variable_info<dim> varInfo;
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


// =====================================================================
// RESIDUAL CONSTRUCTION FUNCTIONS (RHS, LHS, ENERGY DENSITY)
// =====================================================================

template <int dim>
void generalizedProblem<dim>::getRHS(const MatrixFree<dim,double> &data,
					       std::vector<vectorType*> &dst,
					       const std::vector<vectorType*> &src,
					       const std::pair<unsigned int,unsigned int> &cell_range) const{


  //initialize FEEvaulation objects
  std::vector<typeScalar> scalar_vars;
  std::vector<typeVector> vector_vars;

  for (unsigned int i=0; i<num_var; i++){
	  if (varInfoListRHS[i].is_scalar){
		  typeScalar var(data, i);
		  scalar_vars.push_back(var);
	  }
	  else {
		  typeVector var(data, i);
		  vector_vars.push_back(var);
	  }
  }

  std::vector<modelVariable<dim> > modelVarList;
  std::vector<modelResidual<dim> > modelResidualsList;
  modelVarList.reserve(num_var);
  modelResidualsList.reserve(num_var);

  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){

	  // Initialize, read DOFs, and set evaulation flags for each variable
	  for (unsigned int i=0; i<num_var; i++){
		  if (varInfoListRHS[i].is_scalar) {
			  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].reinit(cell);
			  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].read_dof_values_plain(*src[varInfoListRHS[i].global_var_index]);
			  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].evaluate(need_value[i], need_gradient[i], need_hessian[i]);
		  }
		  else {
			  vector_vars[varInfoListRHS[i].scalar_or_vector_index].reinit(cell);
			  vector_vars[varInfoListRHS[i].scalar_or_vector_index].read_dof_values_plain(*src[varInfoListRHS[i].global_var_index]);
			  vector_vars[varInfoListRHS[i].scalar_or_vector_index].evaluate(need_value[i], need_gradient[i], need_hessian[i]);
		  }
	  }

	  unsigned int num_q_points;
	  if (scalar_vars.size() > 0){
		  num_q_points = scalar_vars[0].n_q_points;
	  }
	  else {
		  num_q_points = vector_vars[0].n_q_points;
	  }

	  //loop over quadrature points
	  for (unsigned int q=0; q<num_q_points; ++q){

		  dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc;
		  if (scalar_vars.size() > 0){
			  q_point_loc = scalar_vars[0].quadrature_point(q);
		  }
		  else {
			  q_point_loc = vector_vars[0].quadrature_point(q);
		  }

		  for (unsigned int i=0; i<num_var; i++){
			  if (varInfoListRHS[i].is_scalar) {
				  if (need_value[i]){
					  modelVarList[i].scalarValue = scalar_vars[varInfoListRHS[i].scalar_or_vector_index].get_value(q);
				  }
				  if (need_gradient[i]){
					  modelVarList[i].scalarGrad = scalar_vars[varInfoListRHS[i].scalar_or_vector_index].get_gradient(q);
				  }
				  if (need_hessian[i]){
					  modelVarList[i].scalarHess = scalar_vars[varInfoListRHS[i].scalar_or_vector_index].get_hessian(q);
				  }
			  }
			  else {
				  if (need_value[i]){
					  modelVarList[i].vectorValue = vector_vars[varInfoListRHS[i].scalar_or_vector_index].get_value(q);
				  }
				  if (need_gradient[i]){
					  modelVarList[i].vectorGrad = vector_vars[varInfoListRHS[i].scalar_or_vector_index].get_gradient(q);
				  }
				  if (need_hessian[i]){
					  modelVarList[i].vectorHess = vector_vars[varInfoListRHS[i].scalar_or_vector_index].get_hessian(q);
				  }
			  }
		  }

		  // Calculate the residuals
		  residualRHS(modelVarList,modelResidualsList,q_point_loc);

		  // Submit values
		  for (unsigned int i=0; i<num_var; i++){
			  if (varInfoListRHS[i].is_scalar) {
				  if (value_residual[i] == true){
					  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].submit_value(modelResidualsList[i].scalarValueResidual,q);
				  }
      			  if (gradient_residual[i] == true){
      				  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].submit_gradient(modelResidualsList[i].scalarGradResidual,q);
      			  }
      		  }
      		  else {
      			  if (value_residual[i] == true){
      				  vector_vars[varInfoListRHS[i].scalar_or_vector_index].submit_value(modelResidualsList[i].vectorValueResidual,q);
      			  }
      			  if (gradient_residual[i] == true){
      				  vector_vars[varInfoListRHS[i].scalar_or_vector_index].submit_gradient(modelResidualsList[i].vectorGradResidual,q);
      			  }
      		  }
      	  }

	  }

	  for (unsigned int i=0; i<num_var; i++){
		  if (varInfoListRHS[i].is_scalar) {
			  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].integrate(value_residual[i], gradient_residual[i]);
			  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].distribute_local_to_global(*dst[varInfoListRHS[i].global_var_index]);
		  }
		  else {
			  vector_vars[varInfoListRHS[i].scalar_or_vector_index].integrate(value_residual[i], gradient_residual[i]);
			  vector_vars[varInfoListRHS[i].scalar_or_vector_index].distribute_local_to_global(*dst[varInfoListRHS[i].global_var_index]);
		  }
	  }
  }
}

template <int dim>
void  generalizedProblem<dim>::getLHS(const MatrixFree<dim,double> &data,
					       vectorType &dst,
					       const vectorType &src,
					       const std::pair<unsigned int,unsigned int> &cell_range) const{

	variable_info<dim> resInfoLHS;
	for (unsigned int i=0; i<num_var_LHS; i++){
		if (MatrixFreePDE<dim>::currentFieldIndex == varInfoListLHS[i].global_var_index){
			resInfoLHS = varInfoListLHS[i];
		}
	}

	//initialize FEEvaulation objects
	std::vector<typeScalar> scalar_vars;
	std::vector<typeVector> vector_vars;

	for (unsigned int i=0; i<num_var_LHS; i++){
		if (varInfoListLHS[i].is_scalar){
			typeScalar var(data, varInfoListLHS[i].global_var_index);
			scalar_vars.push_back(var);
		}
		else {
			typeVector var(data, varInfoListLHS[i].global_var_index);
			vector_vars.push_back(var);
		}
	}

	std::vector<modelVariable<dim> > modelVarList;
	modelVarList.reserve(num_var_LHS);
	modelResidual<dim> modelRes;

	//loop over cells
	for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){

		// Initialize, read DOFs, and set evaulation flags for each variable
		for (unsigned int i=0; i<num_var_LHS; i++){
			if (varInfoListLHS[i].is_scalar) {
				scalar_vars[varInfoListLHS[i].scalar_or_vector_index].reinit(cell);
				if ( varInfoListLHS[i].global_var_index == resInfoLHS.global_var_index ){
					scalar_vars[varInfoListLHS[i].scalar_or_vector_index].read_dof_values_plain(src);
				}
				else{
					scalar_vars[varInfoListLHS[i].scalar_or_vector_index].read_dof_values_plain(*MatrixFreePDE<dim>::solutionSet[varInfoListLHS[i].global_var_index]);
				}
				scalar_vars[varInfoListLHS[i].scalar_or_vector_index].evaluate(need_value_LHS[varInfoListLHS[i].global_var_index], need_gradient_LHS[varInfoListLHS[i].global_var_index], need_hessian_LHS[varInfoListLHS[i].global_var_index]);
			}
			else {
				vector_vars[varInfoListLHS[i].scalar_or_vector_index].reinit(cell);
				if ( varInfoListLHS[i].global_var_index == resInfoLHS.global_var_index ){
					vector_vars[varInfoListLHS[i].scalar_or_vector_index].read_dof_values_plain(src);
				}
				else {
					vector_vars[varInfoListLHS[i].scalar_or_vector_index].read_dof_values_plain(*MatrixFreePDE<dim>::solutionSet[varInfoListLHS[i].global_var_index]);
				}
				vector_vars[varInfoListLHS[i].scalar_or_vector_index].evaluate(need_value_LHS[varInfoListLHS[i].global_var_index], need_gradient_LHS[varInfoListLHS[i].global_var_index], need_hessian_LHS[varInfoListLHS[i].global_var_index]);
			}
		}

		unsigned int num_q_points;
		if (scalar_vars.size() > 0){
			num_q_points = scalar_vars[0].n_q_points;
		}
		else {
			num_q_points = vector_vars[0].n_q_points;
		}

		//loop over quadrature points
	    for (unsigned int q=0; q<num_q_points; ++q){
	    	dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc;
	    	if (scalar_vars.size() > 0){
	    		q_point_loc = scalar_vars[0].quadrature_point(q);
	    	}
	    	else {
	    		q_point_loc = vector_vars[0].quadrature_point(q);
	    	}

	    	for (unsigned int i=0; i<num_var_LHS; i++){
	    		if (varInfoListLHS[i].is_scalar) {
	    			if (need_value_LHS[varInfoListLHS[i].global_var_index]){
	    				modelVarList[i].scalarValue = scalar_vars[varInfoListLHS[i].scalar_or_vector_index].get_value(q);
	    			}
	    			if (need_gradient_LHS[varInfoListLHS[i].global_var_index]){
	    				modelVarList[i].scalarGrad = scalar_vars[varInfoListLHS[i].scalar_or_vector_index].get_gradient(q);
	    			}
	    			if (need_hessian_LHS[varInfoListLHS[i].global_var_index]){
	    				modelVarList[i].scalarHess = scalar_vars[varInfoListLHS[i].scalar_or_vector_index].get_hessian(q);
	    			}
	    		}
	    		else {
	    			if (need_value_LHS[varInfoListLHS[i].global_var_index]){
	    				modelVarList[i].vectorValue = vector_vars[varInfoListLHS[i].scalar_or_vector_index].get_value(q);
	    			}
	    			if (need_gradient_LHS[varInfoListLHS[i].global_var_index]){
	    				modelVarList[i].vectorGrad = vector_vars[varInfoListLHS[i].scalar_or_vector_index].get_gradient(q);
	    			}
	    			if (need_hessian_LHS[varInfoListLHS[i].global_var_index]){
	    				modelVarList[i].vectorHess = vector_vars[varInfoListLHS[i].scalar_or_vector_index].get_hessian(q);
	    			}
	    		}
	    	}

	    	// Calculate the residuals
	    	residualLHS(modelVarList,modelRes,q_point_loc);

	    	// Submit values
			if (resInfoLHS.is_scalar){
				if (value_residual[resInfoLHS.global_var_index]){
					scalar_vars[resInfoLHS.scalar_or_vector_index].submit_value(modelRes.scalarValueResidual,q);
				}
				if (gradient_residual[resInfoLHS.global_var_index]){
					scalar_vars[resInfoLHS.scalar_or_vector_index].submit_gradient(modelRes.scalarGradResidual,q);
				}
			}
			else {
				if (value_residual[resInfoLHS.global_var_index]){
					vector_vars[resInfoLHS.scalar_or_vector_index].submit_value(modelRes.vectorValueResidual,q);
				}
				if (gradient_residual[resInfoLHS.global_var_index]){
					vector_vars[resInfoLHS.scalar_or_vector_index].submit_gradient(modelRes.vectorGradResidual,q);
				}
			}

	    }

	    //integrate
		if (resInfoLHS.is_scalar) {
			scalar_vars[resInfoLHS.scalar_or_vector_index].integrate(value_residual[resInfoLHS.global_var_index], gradient_residual[resInfoLHS.global_var_index]);
			scalar_vars[resInfoLHS.scalar_or_vector_index].distribute_local_to_global(dst);
		}
		else {
			vector_vars[resInfoLHS.scalar_or_vector_index].integrate(value_residual[resInfoLHS.global_var_index], gradient_residual[resInfoLHS.global_var_index]);
			vector_vars[resInfoLHS.scalar_or_vector_index].distribute_local_to_global(dst);
		}
	}

}

// Calculate the free energy
template <int dim>
void  generalizedProblem<dim>::getEnergy(const MatrixFree<dim,double> &data,
				    std::vector<vectorType*> &dst,
				    const std::vector<vectorType*> &src,
				    const std::pair<unsigned int,unsigned int> &cell_range) {

	//initialize FEEvaulation objects
	  std::vector<typeScalar> scalar_vars;
	  std::vector<typeVector> vector_vars;

	  for (unsigned int i=0; i<num_var; i++){
		  if (varInfoListRHS[i].is_scalar){
			  typeScalar var(data, i);
			  scalar_vars.push_back(var);
		  }
		  else {
			  typeVector var(data, i);
			  vector_vars.push_back(var);
		  }
	  }

	  std::vector<modelVariable<dim> > modelVarList;
	  std::vector<modelResidual<dim> > modelResidualsList;
	  modelVarList.reserve(num_var);
	  modelResidualsList.reserve(num_var);

	  //loop over cells
	  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){

		  // Initialize, read DOFs, and set evaulation flags for each variable
		  for (unsigned int i=0; i<num_var; i++){
			  if (varInfoListRHS[i].is_scalar) {
				  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].reinit(cell);
				  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].read_dof_values_plain(*src[varInfoListRHS[i].global_var_index]);
				  scalar_vars[varInfoListRHS[i].scalar_or_vector_index].evaluate(need_value[i], need_gradient[i], need_hessian[i]);
			  }
			  else {
				  vector_vars[varInfoListRHS[i].scalar_or_vector_index].reinit(cell);
				  vector_vars[varInfoListRHS[i].scalar_or_vector_index].read_dof_values_plain(*src[varInfoListRHS[i].global_var_index]);
				  vector_vars[varInfoListRHS[i].scalar_or_vector_index].evaluate(need_value[i], need_gradient[i], need_hessian[i]);
			  }
		  }


		  unsigned int num_q_points;
		  if (scalar_vars.size() > 0){
			  num_q_points = scalar_vars[0].n_q_points;
		  }
		  else {
			  num_q_points = vector_vars[0].n_q_points;
		  }

		  dealii::AlignedVector<dealii::VectorizedArray<double> > JxW(num_q_points);

		  if (scalar_vars.size() > 0){
			  scalar_vars[0].fill_JxW_values(JxW);
		  }
		  else {
			  vector_vars[0].fill_JxW_values(JxW);
		  }

		  //loop over quadrature points
		  for (unsigned int q=0; q<num_q_points; ++q){
			  dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc;
			  if (scalar_vars.size() > 0){
				  q_point_loc = scalar_vars[0].quadrature_point(q);
			  }
			  else {
				  q_point_loc = vector_vars[0].quadrature_point(q);
			  }

			  for (unsigned int i=0; i<num_var; i++){
				  if (varInfoListRHS[i].is_scalar) {
					  if (need_value[i]){
						  modelVarList[i].scalarValue = scalar_vars[varInfoListRHS[i].scalar_or_vector_index].get_value(q);
					  }
					  if (need_gradient[i]){
						  modelVarList[i].scalarGrad = scalar_vars[varInfoListRHS[i].scalar_or_vector_index].get_gradient(q);
					  }
					  if (need_hessian[i]){
						  modelVarList[i].scalarHess = scalar_vars[varInfoListRHS[i].scalar_or_vector_index].get_hessian(q);
					  }
				  }
				  else {
					  if (need_value[i]){
						  modelVarList[i].vectorValue = vector_vars[varInfoListRHS[i].scalar_or_vector_index].get_value(q);
					  }
					  if (need_gradient[i]){
						  modelVarList[i].vectorGrad = vector_vars[varInfoListRHS[i].scalar_or_vector_index].get_gradient(q);
					  }
					  if (need_hessian[i]){
						  modelVarList[i].vectorHess = vector_vars[varInfoListRHS[i].scalar_or_vector_index].get_hessian(q);
					  }
				  }
			  }

			  // Calculate the energy density
			  energyDensity(modelVarList,JxW[q],q_point_loc);
		  }
	  }

}

//compute the integral of one of the fields
template <int dim>
void generalizedProblem<dim>::computeIntegral(double& integratedField){
  QGauss<dim>  quadrature_formula(finiteElementDegree+1);
  FE_Q<dim> FE (QGaussLobatto<1>(finiteElementDegree+1));
  FEValues<dim> fe_values (FE, quadrature_formula, update_values | update_JxW_values | update_quadrature_points);
  const unsigned int   dofs_per_cell = FE.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();
  std::vector<double> cVal(n_q_points);

  typename DoFHandler<dim>::active_cell_iterator cell= this->dofHandlersSet[0]->begin_active(), endc = this->dofHandlersSet[0]->end();

  double value = 0.0;

  unsigned int fieldIndex;
  fieldIndex=this->getFieldIndex("c");

  for (; cell!=endc; ++cell) {
	  if (cell->is_locally_owned()){
    	fe_values.reinit (cell);

    	fe_values.get_function_values(*this->solutionSet[fieldIndex], cVal);

    	for (unsigned int q=0; q<n_q_points; ++q){
    		value+=(cVal[q])*fe_values.JxW(q);
    	}
	  }
  }

  value=Utilities::MPI::sum(value, MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0){
  std::cout<<"Integrated field: "<<value<<std::endl;
  }

  integratedField = value;
}

// =====================================================================
// ADAPTIVE MESHING FUNCTIONS
// =====================================================================

//adaptive refinement control
template <int dim>
void generalizedProblem<dim>::adaptiveRefine(unsigned int currentIncrement){
	#if hAdaptivity == true
	//if ((currentIncrement>0) && (currentIncrement%skipRemeshingSteps==0)){
//	if ( (currentIncrement == 2) || (currentIncrement%skipRemeshingSteps==0) ){
//		this->refineMesh(currentIncrement);
//	}

	if ( (currentIncrement == 0) ){
		for (unsigned int remesh_index=0; remesh_index < (maxRefinementLevel-minRefinementLevel); remesh_index++){
			this->refineMesh(currentIncrement);
		}
	}
	else if ( (currentIncrement%skipRemeshingSteps==0) ){
			this->refineMesh(currentIncrement);
		}
	#endif
}

//adaptive refinement criterion
template <int dim>
void generalizedProblem<dim>::adaptiveRefineCriterion(){
#if hAdaptivity == true
	//Custom defined estimation criterion
	std::vector<int> refine_criterion_fields = refineCriterionFields;
	std::vector<double> refine_window_max = refineWindowMax;
	std::vector<double> refine_window_min = refineWindowMin;
	std::vector<std::vector<double> > errorOutV;


	QGauss<dim>  quadrature(finiteElementDegree+1);
	FEValues<dim> fe_values (*this->FESet[refine_criterion_fields[0]], quadrature, update_values);
	const unsigned int   num_quad_points = quadrature.size();

	std::vector<double> errorOut(num_quad_points);

	typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandlersSet_nonconst[refine_criterion_fields[0]]->begin_active(), endc = this->dofHandlersSet_nonconst[refine_criterion_fields[0]]->end();

	for (;cell!=endc; ++cell){
		if (cell->is_locally_owned()){
			fe_values.reinit (cell);

			for (unsigned int field_index=0; field_index<refine_criterion_fields.size(); field_index++){
				fe_values.get_function_values(*this->solutionSet[refine_criterion_fields[field_index]], errorOut);
				errorOutV.push_back(errorOut);
			}

			bool mark_refine = false;

			for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
				for (unsigned int field_index=0; field_index<refine_criterion_fields.size(); field_index++){
					if ((errorOutV[field_index][q_point]>refine_window_min[field_index]) && (errorOutV[field_index][q_point]<refine_window_max[field_index])){
						mark_refine = true;
						break;
					}
				}
			}

			errorOutV.clear();

//			fe_values.get_function_values(*this->solutionSet[refine_criterion_fields[0]], errorOut);
//
//			bool mark_refine = false;
//
//			for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
//				if ((errorOut[q_point]>refine_window_min[0]) && (errorOut[q_point]<refine_window_max[0])){
//					mark_refine = true;
//					break;
//				}
//			}


			if ( (mark_refine == true) ){
				cell->set_refine_flag();
			}
			else {
				cell->set_coarsen_flag();
			}
		}
	}
#endif
}


// =====================================================================
// FUNCTION TO BUILD THE VECTOR OF FIELDS
// =====================================================================

template <int dim>
void generalizedProblem<dim>::buildFields(){
	// Build each of the fields in the system
	for (unsigned int i=0; i<num_var; i++){
		  if (var_type[i] == "SCALAR"){
			  if (var_eq_type[i] == "ELLIPTIC"){
				  this->fields.push_back(Field<problemDIM>(SCALAR, ELLIPTIC, var_name[i]));
			  }
			  else if (var_eq_type[i] == "PARABOLIC"){
				  this->fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, var_name[i]));
			  }
			  else{
				  // Need to change to throw an exception
				  std::cerr << "Error: Equation type must be ELLIPTIC or PARABOLIC " << std::endl;
			  }
		  }
		  else if (var_type[i] == "VECTOR"){
			  if (var_eq_type[i] == "ELLIPTIC"){
				  this->fields.push_back(Field<problemDIM>(VECTOR, ELLIPTIC, var_name[i]));
			  }
			  else if (var_eq_type[i] == "PARABOLIC"){
				  this->fields.push_back(Field<problemDIM>(VECTOR, PARABOLIC, var_name[i]));
			  }
			  else{
				  // Need to change to throw an exception
				  std::cerr << "Error: Variable type must be SCALAR or VECTOR " << std::endl;
			  }
		  }
	  }

}

// =====================================================================
// INITIAL CONDITION FUNCTIONS
// =====================================================================

#ifdef enablePFields
#if enablePFields == true
template <int dim>
class InitialConditionPField : public Function<dim>
{
public:
  unsigned int index;
  Vector<double> values;
  typedef PRISMS::PField<double*, double, 2> ScalarField2D;
  ScalarField2D &inputField;

  InitialConditionPField (const unsigned int _index, ScalarField2D &_inputField) : Function<dim>(1), index(_index), inputField(_inputField) {}

  double value (const Point<dim> &p, const unsigned int component = 0) const
  {
	  double scalar_IC;

	  double coord[dim];
	  for (unsigned int i = 0; i < dim; i++){
		  coord[i] = p(i);
	  }

	  scalar_IC = inputField(coord);

	  return scalar_IC;
  }
};

#endif
#else
#define enablePFields false
#endif

// =================================================================================

//apply initial conditions
template <int dim>
void generalizedProblem<dim>::applyInitialConditions()
{


#ifndef loadICs
#define loadICs {}
#endif

#ifndef loadSerialFile
#define loadSerialFile {}
#endif

#ifndef loadFileName
#define loadFileName {}
#endif

#ifndef loadFieldName
#define loadFieldName {}
#endif

std::vector<bool> load_ICs = loadICs;
std::vector<bool> load_serial_file = loadSerialFile;
std::vector<std::string> load_file_name = loadFileName;
std::vector<std::string> load_field_name = loadFieldName;

// If load_ICs is empty, it should be set to false for each variable
if (load_ICs.size() == 0){
	for (unsigned int i=0; i<num_var; i++){
		load_ICs.push_back(false);
	}
}

for (unsigned int var_index=0; var_index < num_var; var_index++){
	if (load_ICs[var_index] == false){
		if (var_type[var_index] == "SCALAR"){
			VectorTools::interpolate (*this->dofHandlersSet[var_index], InitialCondition<dim>(var_index), *this->solutionSet[var_index]);
		}
		else {
			VectorTools::interpolate (*this->dofHandlersSet[var_index], InitialConditionVec<dim>(var_index), *this->solutionSet[var_index]);
		}
	}
	else{
		#if enablePFields == true
		// Declare the PField types and containers
		typedef PRISMS::PField<double*, double, 2> ScalarField2D;
		typedef PRISMS::Body<double*, 2> Body2D;
		Body2D body;

		// Create the filename of the the file to be loaded
		std::string filename;
		if (load_serial_file[var_index] == true){
			filename = load_file_name[var_index] + ".vtk";
		}
		else {
			int proc_num = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
			std::ostringstream conversion;
			conversion << proc_num;
			filename = load_file_name[var_index] + "." + conversion.str() + ".vtk";
		}

		// Load the data from the file using a PField
		body.read_vtk(filename);
		ScalarField2D &conc = body.find_scalar_field(load_field_name[var_index]);
		if (var_type[var_index] == "SCALAR"){
			VectorTools::interpolate (*this->dofHandlersSet[var_index], InitialConditionPField<dim>(var_index,conc), *this->solutionSet[var_index]);
		}
		else {
			std::cout << "PRISMS-PF Error: Cannot load vector fields. Loading initial conditions from file is currently limited to scalar fields" << std::endl;
		}
		#else
		std::cout << "PRISMS-PF Error: The parameter \"enablePFields\" must be set to true to load initial conditions from file." << std::endl;
		#endif
	}
}
}

// =================================================================================

// I don't think vector fields are implemented in PFields yet
//template <int dim>
//class InitialConditionPFieldVec : public Function<dim>
//{
//public:
//  unsigned int index;
//  Vector<double> values;
//  typedef PRISMS::PField<double*, double, 2> ScalarField2D;
//  ScalarField2D &inputField;
//
//  InitialConditionPFieldVec (const unsigned int _index, ScalarField2D &_inputField) : Function<dim>(1), index(_index), inputField(_inputField) {}
//
//  void vector_value (const Point<dim> &p,Vector<double> &vector_IC) const
//  {
//	  double coord[dim];
//	  for (unsigned int i = 0; i < dim; i++){
//		  coord[i] = p(i);
//	  }
//
//	  vector_IC = inputField(coord);
//  }
//};


// =====================================================================
// BOUNDARY CONDITION FUNCTIONS
// =====================================================================

template <int dim>
class vectorBCFunction : public Function<dim,double>
{
public:
	vectorBCFunction(const std::vector<double> BC_values);
  virtual void vector_value (const Point<dim> &p, Vector<double>   &values) const;

  virtual void vector_value_list (const std::vector<Point<dim> > &points, std::vector<Vector<double> > &value_list) const;

private:
  const std::vector<double> BC_values;
};

template <int dim>
vectorBCFunction<dim>::vectorBCFunction(std::vector<double> input_values) : Function<dim>(dim), BC_values (input_values) {}

template <int dim>
void vectorBCFunction<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const {

	for (unsigned int i=0; i<dim; i++) {
		values(i) = BC_values[i];
	}
}

template <int dim>
void vectorBCFunction<dim>::vector_value_list (const std::vector<Point<dim> > &points, std::vector<Vector<double> > &value_list) const{
	const unsigned int n_points = points.size();
	for (unsigned int p=0; p<n_points; ++p)
		vectorBCFunction<dim>::vector_value(points[p],value_list[p]);
}

//apply Dirchlet BC function
template <int dim>
void generalizedProblem<dim>::applyDirichletBCs(){

  // First, get the variable index of the current field
  unsigned int starting_BC_list_index = 0;

  for (unsigned int i=0; i<this->currentFieldIndex; i++){

	  if (var_type[i] == "SCALAR"){
		  starting_BC_list_index++;
	  }
	  else {
		  starting_BC_list_index+=dim;
	  }
  }

  if (var_type[this->currentFieldIndex] == "SCALAR"){
	  for (unsigned int direction = 0; direction < 2*dim; direction++){
		  if (BC_list[starting_BC_list_index].var_BC_type[direction] == "DIRICHLET"){
			  VectorTools::interpolate_boundary_values (*this->dofHandlersSet[this->currentFieldIndex],\
					  direction, ConstantFunction<dim>(BC_list[starting_BC_list_index].var_BC_val[direction],1), *(ConstraintMatrix*) \
					  this->constraintsDirichletSet[this->currentFieldIndex]);
		  }
	  }
  }
  else {
	  for (unsigned int direction = 0; direction < 2*dim; direction++){

		  std::vector<double> BC_values;
		  for (unsigned int component=0; component < dim; component++){
			  BC_values.push_back(BC_list[starting_BC_list_index+component].var_BC_val[direction]);
		  }

		  std::vector<bool> mask;
		  for (unsigned int component=0; component < dim; component++){
			  if (BC_list[starting_BC_list_index+component].var_BC_type[direction] == "DIRICHLET"){
				  mask.push_back(true);
			  }
			  else {
				  mask.push_back(false);
			  }
		  }

		  VectorTools::interpolate_boundary_values (*this->dofHandlersSet[this->currentFieldIndex],\
				  direction, vectorBCFunction<dim>(BC_values), *(ConstraintMatrix*) \
				  this->constraintsDirichletSet[this->currentFieldIndex],mask);


	  }
  }
}

//methods to mark boundaries
template <int dim>
void generalizedProblem<dim>::markBoundaries(){

	std::vector<double> domain_size;
	domain_size.push_back(spanX);
	domain_size.push_back(spanY);
	domain_size.push_back(spanZ);

	typename Triangulation<dim>::cell_iterator
	cell = MatrixFreePDE<dim>::triangulation.begin (),
	endc = MatrixFreePDE<dim>::triangulation.end();

	for (; cell!=endc; ++cell){

		// Mark all of the faces
		for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell;++face_number){
			for (unsigned int i=0; i<dim; i++){
				if ( std::fabs(cell->face(face_number)->center()(i) - (0)) < 1e-12 ){
					cell->face(face_number)->set_boundary_id (2*i);
				}
				else if (std::fabs(cell->face(face_number)->center()(i) - (domain_size[i])) < 1e-12){
					cell->face(face_number)->set_boundary_id (2*i+1);
				}

			}
		}
	}
}

// Input the boundary conditions for each face individually for 3D domains
template <int dim>
void generalizedProblem<dim>::inputBCs(int var, int component, std::string BC_type_dim1_min, double BC_value_dim1_min,
		std::string BC_type_dim1_max, double BC_value_dim1_max, std::string BC_type_dim2_min, double BC_value_dim2_min,
		std::string BC_type_dim2_max, double BC_value_dim2_max,std::string BC_type_dim3_min, double BC_value_dim3_min,
		std::string BC_type_dim3_max, double BC_value_dim3_max){

	// Validate input
	try{
		if ((BC_type_dim1_min == "PERIODIC") && (BC_type_dim1_max != "PERIODIC")){
			throw 0;
		}
		if ((BC_type_dim2_min == "PERIODIC") && (BC_type_dim2_max != "PERIODIC")){
			throw 0;
		}
		if ((BC_type_dim3_min == "PERIODIC") && (BC_type_dim3_max != "PERIODIC")){
			throw 0;
		}
	}
	catch (int e){
		if (e == 0){
			std::cout << "Error: For periodic BCs, both faces for a given direction must be set as periodic. "
					"Please check the BCs that are set in ICs_and_BCs.h." << std::endl;
		}
		abort();
	}

	varBCs<dim> newBC;
	newBC.var_BC_type.push_back(BC_type_dim1_min);
	newBC.var_BC_type.push_back(BC_type_dim1_max);
	newBC.var_BC_type.push_back(BC_type_dim2_min);
	newBC.var_BC_type.push_back(BC_type_dim2_max);
	newBC.var_BC_type.push_back(BC_type_dim3_min);
	newBC.var_BC_type.push_back(BC_type_dim3_max);

	newBC.var_BC_val.push_back(BC_value_dim1_min);
	newBC.var_BC_val.push_back(BC_value_dim1_max);
	newBC.var_BC_val.push_back(BC_value_dim2_min);
	newBC.var_BC_val.push_back(BC_value_dim2_max);
	newBC.var_BC_val.push_back(BC_value_dim3_min);
	newBC.var_BC_val.push_back(BC_value_dim3_max);

	BC_list.push_back(newBC);
}

// Input the boundary conditions for each face individually for 2D domains
template <int dim>
void generalizedProblem<dim>::inputBCs(int var, int component, std::string BC_type_dim1_min, double BC_value_dim1_min,
		std::string BC_type_dim1_max, double BC_value_dim1_max, std::string BC_type_dim2_min, double BC_value_dim2_min,
		std::string BC_type_dim2_max, double BC_value_dim2_max){

	// Validate input
	try{
		if ((BC_type_dim1_min == "PERIODIC") && (BC_type_dim1_max != "PERIODIC")){
			throw 0;
		}
		if ((BC_type_dim2_min == "PERIODIC") && (BC_type_dim2_max != "PERIODIC")){
			throw 0;
		}
	}
	catch (int e){
		if (e == 0){
			std::cout << "Error: For periodic BCs, both faces for a given direction must be set as periodic. "
					"Please check the BCs that are set in ICs_and_BCs.h." << std::endl;
		}
		abort();
	}

	varBCs<dim> newBC;
	newBC.var_BC_type.push_back(BC_type_dim1_min);
	newBC.var_BC_type.push_back(BC_type_dim1_max);
	newBC.var_BC_type.push_back(BC_type_dim2_min);
	newBC.var_BC_type.push_back(BC_type_dim2_max);

	newBC.var_BC_val.push_back(BC_value_dim1_min);
	newBC.var_BC_val.push_back(BC_value_dim1_max);
	newBC.var_BC_val.push_back(BC_value_dim2_min);
	newBC.var_BC_val.push_back(BC_value_dim2_max);

	BC_list.push_back(newBC);
}

// Input the boundary conditions for each face individually for 1D domains
template <int dim>
void generalizedProblem<dim>::inputBCs(int var, int component, std::string BC_type_dim1_min, double BC_value_dim1_min,
		std::string BC_type_dim1_max, double BC_value_dim1_max){

	// Validate input
	try{
		if ((BC_type_dim1_min == "PERIODIC") && (BC_type_dim1_max != "PERIODIC")){
			throw 0;
		}

	}
	catch (int e){
		if (e == 0){
			std::cout << "Error: For periodic BCs, both faces for a given direction must be set as periodic. "
					"Please check the BCs that are set in ICs_and_BCs.h." << std::endl;
		}
		abort();
	}


	varBCs<dim> newBC;
	newBC.var_BC_type.push_back(BC_type_dim1_min);
	newBC.var_BC_type.push_back(BC_type_dim1_max);

	newBC.var_BC_val.push_back(BC_value_dim1_min);
	newBC.var_BC_val.push_back(BC_value_dim1_max);

	BC_list.push_back(newBC);
}

// Input the boundary conditions when all faces have the same boundary condition
template <int dim>
void generalizedProblem<dim>::inputBCs(int var, int component, std::string BC_type, double BC_value){

	varBCs<dim> newBC;
	for (unsigned int face=0; face<(dim*2); face++){
		newBC.var_BC_type.push_back(BC_type);
		newBC.var_BC_val.push_back(BC_value);
	}

	BC_list.push_back(newBC);
}

// Based on the contents of BC_list, mark faces on the triangulation as periodic
template <int dim>
void generalizedProblem<dim>::setPeriodicity(){
	std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> > periodicity_vector;
	for (int i=0; i<dim; ++i){
		bool periodic_pair = false;
		for (unsigned int field_num=0; field_num < BC_list.size(); field_num++){
			if (BC_list[field_num].var_BC_type[2*i] == "PERIODIC"){
				periodic_pair = true;
			}
		}
		if (periodic_pair == true){
			GridTools::collect_periodic_faces(this->triangulation, /*b_id1*/ 2*i, /*b_id2*/ 2*i+1,
							/*direction*/ i, periodicity_vector);
		}
	}

	this->triangulation.add_periodicity(periodicity_vector);
	this->pcout << "periodic facepairs: " << periodicity_vector.size() << std::endl;
}

// Set constraints to enforce periodic boundary conditions
template <int dim>
void generalizedProblem<dim>::setPeriodicityConstraints(ConstraintMatrix * constraints, DoFHandler<dim>* dof_handler){
	// First, get the variable index of the current field
	unsigned int starting_BC_list_index = 0;
	for (unsigned int i=0; i<this->currentFieldIndex; i++){
		if (var_type[i] == "SCALAR"){
			starting_BC_list_index++;
		}
		else {
			starting_BC_list_index+=dim;
		}
	}

	std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator> > periodicity_vector;
    for (int i=0; i<dim; ++i){
    	if (BC_list[starting_BC_list_index].var_BC_type[2*i] == "PERIODIC"){
    		GridTools::collect_periodic_faces(*dof_handler, /*b_id1*/ 2*i, /*b_id2*/ 2*i+1,
    				/*direction*/ i, periodicity_vector);
    	}
    }
    DoFTools::make_periodicity_constraints<DoFHandler<dim> >(periodicity_vector, *constraints);
}

// Determine which (if any) components of the current field have rigid body modes (i.e no Dirichlet BCs) if the
// equation is elliptic
template <int dim>
void generalizedProblem<dim>::getComponentsWithRigidBodyModes(std::vector<int> & rigidBodyModeComponents){

	// Rigid body modes only matter for elliptic equations
	if (var_eq_type[this->currentFieldIndex] == "ELLIPTIC"){

		// First, get the variable index of the current field
		unsigned int starting_BC_list_index = 0;
		for (unsigned int i=0; i<this->currentFieldIndex; i++){
			if (var_type[i] == "SCALAR"){
				starting_BC_list_index++;
			}
			else {
				starting_BC_list_index+=dim;
			}
		}

		// Get number of components of the field
		unsigned int num_components = 1;
		if (var_type[this->currentFieldIndex] == "VECTOR"){
			num_components = dim;
		}

		// Loop over each component and determine if it has a rigid body mode (i.e. no Dirichlet BCs)
		for (unsigned int component=0; component < num_components; component++){
			bool rigidBodyMode = true;
			for (unsigned int direction = 0; direction < 2*dim; direction++){

				if (BC_list[starting_BC_list_index+component].var_BC_type[direction] == "DIRICHLET"){
					rigidBodyMode = false;
				}

			}
			// If the component has a rigid body mode, add it to the list
			if (rigidBodyMode == true){
				rigidBodyModeComponents.push_back(component);
			}
		}
	}
}

// =====================================================================
// NUCLEATION FUNCTIONS
// =====================================================================

//structure representing each nucleus
struct nucleus{
  unsigned int index;
  dealii::Point<problemDIM> center;
  double radius;
  double seededTime, seedingTime;
};
//vector of all nucleus seeded in the problem
std::vector<nucleus> nuclei, localNuclei;

//nucleation model implementation
template <int dim>
void generalizedProblem<dim>::modifySolutionFields()
{
  //current time
  double t=this->currentTime;
  unsigned int inc=this->currentIncrement;
  double dx=spanX/( (double)subdivisionsX )/std::pow(2.0,refineFactor);
  double rand_val;
  int count = 0;
  //nucleation parameters
  double nRadius = 2.5; //spanX/20.0;
  double minDistBetwenNuclei=4*nRadius;

  unsigned int maxNumberNuclei=5; // doesn't do anything currently

  //get the list of node points in the domain
  std::map<dealii::types::global_dof_index, dealii::Point<dim> > support_points;
  dealii::DoFTools::map_dofs_to_support_points (dealii::MappingQ1<dim>(), *this->dofHandlersSet[0], support_points);
  //fields
  vectorType* n1=this->solutionSet[this->getFieldIndex("n1")];
  vectorType* n2=this->solutionSet[this->getFieldIndex("n2")];
  vectorType* n3=this->solutionSet[this->getFieldIndex("n3")];
  vectorType* c=this->solutionSet[this->getFieldIndex("c")];
  const double k1 = 0.0001; // nucleation probability constant
  const double k2 = 1.0;	// nucleation probability constant
  const double c0 = 0.300;	// baseline concentration?
  double J = 0.0;
  //delete the previous entries in the nuclei vector (old nucleus are still retained in the localNuclei vector)
  nuclei.clear();

#ifndef c_matrix
  double c_matrix = 1.0e-6;
#endif

  //populate localNuclei vector
  if (inc <= timeIncrements){
    nucleus* temp;
    //add nuclei based on concentration field values
    //loop over all points in the domain
    for (typename std::map<dealii::types::global_dof_index, dealii::Point<dim> >::iterator it=support_points.begin(); it!=support_points.end(); ++it){
      unsigned int dof=it->first;
      //set only local owned values of the parallel vector (eventually turn this into a separate function for each order parameter)
      if (n1->locally_owned_elements().is_element(dof)){
    	  dealii::Point<dim> nodePoint=it->second;
    	  double n1Value=(*n1)(dof);
    	  double n2Value=(*n2)(dof);
    	  double n3Value=(*n3)(dof);
    	  double cValue=(*c)(dof);

    	  rand_val = (double)rand()/(RAND_MAX);

    	  if ((t > 1000000000*timeStep) || (n1Value+n2Value+n3Value > 1.0e-6) || (cValue <= 0.0)) {
    		  J = 0;
    	  }
		  else{
			  J = cValue/c_matrix * dx*dx/((double)spanX * (double)spanY) * 0.01; // Only true in 2D!
    	  }

    	  if (rand_val <= J){
    		  bool isClose=false;
    		  for (std::vector<nucleus>::iterator thisNuclei=localNuclei.begin(); thisNuclei!=localNuclei.end(); ++thisNuclei){
    			  if (thisNuclei->center.distance(nodePoint)<minDistBetwenNuclei){
    				  isClose=true;
    			  }
    		  }

    		  if (!isClose){
    			  temp = new nucleus;
    			  temp->index=localNuclei.size();
    			  temp->center=nodePoint;
    			  temp->radius=nRadius;
    			  temp->seededTime=t;
    			  temp->seedingTime = 10000.0*timeStep;
    			  localNuclei.push_back(*temp);
    		  }
    	  }
      }
    }


    //filter nuclei by comparing with other processors
    int numProcs=Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
    int thisProc=Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
    std::vector<int> numNucleiInProcs(numProcs, 0);
    //send nuclei information to processor 0
    int numNuclei=localNuclei.size();
    //send information about number of nuclei to processor 0
    if (thisProc!=0){
      MPI_Send(&numNuclei, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }
    else{
      numNucleiInProcs[0]=numNuclei;
      for (int proc=1; proc<numProcs; proc++){
	MPI_Recv(&numNucleiInProcs[proc], 1, MPI_INT, proc, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //filter nuclei in processor zero
    //receive nuclei info from all processors
    if (thisProc!=0){
    	if (numNuclei>0){
    		std::vector<double> tempData((dim+3)*numNuclei);
    		unsigned int i=0;
    		for (std::vector<nucleus>::iterator thisNuclei=localNuclei.begin(); thisNuclei!=localNuclei.end(); ++thisNuclei){
    			tempData[i*(dim+3)]=thisNuclei->radius;
    			tempData[i*(dim+3)+1]=thisNuclei->seededTime;
    			tempData[i*(dim+3)+2]=thisNuclei->seedingTime;
    			for (unsigned int j=0; j<dim; j++) tempData[i*(dim+3)+3+j]=thisNuclei->center[j];
    			i++;
    		}
    		MPI_Send(&tempData[0], numNuclei*(dim+3), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
    	}
    }
    else{
    	//temporary array to store all the nuclei
    	std::vector<std::vector<double>*> tempNuceli(numProcs);
    	for (int proc=0; proc<numProcs; proc++) {
    		std::vector<double>* temp=new std::vector<double>(numNucleiInProcs[proc]*(dim+3));
    		if (numNucleiInProcs[proc]>0){
    			if (proc==0){
    				unsigned int i=0;
    				for (std::vector<nucleus>::iterator thisNuclei=localNuclei.begin(); thisNuclei!=localNuclei.end(); ++thisNuclei){
    					(*temp)[i*(dim+3)]=thisNuclei->radius;
    					(*temp)[i*(dim+3)+1]=thisNuclei->seededTime;
    					(*temp)[i*(dim+3)+2]=thisNuclei->seedingTime;
    					for (unsigned int j=0; j<dim; j++) (*temp)[i*(dim+3)+3+j]=thisNuclei->center[j];
    					i++;
    				}
    			}
    			else{
    				MPI_Recv(&((*temp)[0]), numNucleiInProcs[proc]*(dim+3), MPI_DOUBLE, proc, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    			}
    			tempNuceli[proc]=temp;
    		}
    	}

    	//filter the nuclei and add to nuclei vector in processor zero
    	for (int proc1=0; proc1<numProcs; proc1++) {
    		for (int i1=0; i1<numNucleiInProcs[proc1]; i1++){
    			double rad1=(*tempNuceli[proc1])[i1*(dim+3)];
    			double time1=(*tempNuceli[proc1])[i1*(dim+3)+1];
    			double seedingTime1=(*tempNuceli[proc1])[i1*(dim+3)+2];
    			dealii::Point<dim> center1;
    			for (unsigned int j1=0; j1<dim; j1++) {
    				center1[j1]=(*tempNuceli[proc1])[i1*(dim+3)+3+j1];
    			}
    			bool addNuclei=true;
    			//check if this nuceli present in any other processor
    			for (int proc2=0; proc2<numProcs; proc2++) {
    				if (proc1!=proc2){
    					for (int i2=0; i2<numNucleiInProcs[proc2]; i2++){
    						double rad2=(*tempNuceli[proc2])[i2*(dim+3)];
    						double time2=(*tempNuceli[proc2])[i2*(dim+3)+1];
    						dealii::Point<dim> center2;
    						for (unsigned int j2=0; j2<dim; j2++) center2(j2)=(*tempNuceli[proc2])[i2*(dim+3)+3+j2];
    						if ((center1.distance(center2)<=minDistBetwenNuclei) && (time1>=time2)){
    							addNuclei=false;
    							break;
    						}
    					}
    					if (!addNuclei) {break;}
    				}
    			}
    			if (addNuclei){
    				temp = new nucleus;
    				temp->index=nuclei.size();
    				temp->radius=rad1;
    				temp->seededTime=time1;
    				temp->seedingTime=seedingTime1;
    				temp->center=center1;
    				nuclei.push_back(*temp);
    			}
    		}
    	}
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //disperse nuclei to all other processors
    unsigned int numGlobalNuclei;
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0) {numGlobalNuclei=nuclei.size();}
    MPI_Bcast(&numGlobalNuclei, 1, MPI_INT, 0, MPI_COMM_WORLD);
    this->pcout << "total number of nuclei currently seeded : "  << numGlobalNuclei << std::endl;
    MPI_Barrier(MPI_COMM_WORLD);
    //
    std::vector<double> temp2(numGlobalNuclei*(dim+3));
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0){
      unsigned int i=0;
      for (std::vector<nucleus>::iterator thisNuclei=nuclei.begin(); thisNuclei!=nuclei.end(); ++thisNuclei){
	temp2[i*(dim+3)]=thisNuclei->radius;
	temp2[i*(dim+3)+1]=thisNuclei->seededTime;
	temp2[i*(dim+3)+2]=thisNuclei->seedingTime;
	for (unsigned int j=0; j<dim; j++) temp2[i*(dim+3)+3+j]=thisNuclei->center[j];
	i++;
      }
    }
    MPI_Bcast(&temp2[0], numGlobalNuclei*(dim+3), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    //receive all nuclei
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)!=0){
    	for(unsigned int i=0; i<numGlobalNuclei; i++){
    		temp = new nucleus;
    		temp->index=nuclei.size();
    		temp->radius=temp2[i*(dim+3)];
    		temp->seededTime=temp2[i*(dim+3)+1];
    		temp->seedingTime=temp2[i*(dim+3)+2];
    		dealii::Point<dim> tempCenter;
    		for (unsigned int j=0; j<dim; j++) tempCenter[j]=temp2[i*(dim+3)+3+j];
    		temp->center=tempCenter;
    		nuclei.push_back(*temp);
      }
    }
  }

  //seed nuclei
  unsigned int fieldIndex=this->getFieldIndex("n1");
  for (std::vector<nucleus>::iterator thisNuclei=nuclei.begin(); thisNuclei!=nuclei.end(); ++thisNuclei){

	  dealii::Point<dim> center=thisNuclei->center;
	  double radius=thisNuclei->radius;
	  double seededTime=thisNuclei->seededTime;
	  double seedingTime=thisNuclei->seedingTime;
	  this->pcout << "times: " << t << " " << seededTime << " " << seedingTime << std::endl;
	  //loop over all points in the domain
	  for (typename std::map<dealii::types::global_dof_index, dealii::Point<dim> >::iterator it=support_points.begin(); it!=support_points.end(); ++it){
		  unsigned int dof=it->first;
		  //set only local owned values of the parallel vector
		  if (n1->locally_owned_elements().is_element(dof)){
			  dealii::Point<dim> nodePoint=it->second;
			  //check conditions and seed nuclei
			  double r=nodePoint.distance(center);
			  if (r<=(2*radius)){
				  if ((t>seededTime) && (t<(seededTime+seedingTime))){
					  //this->pcout << "times: " << t << " " << seededTime << " " << seedingTime << std::endl;
					  //(*n1)(dof)=0.5*(1.0-std::tanh((r-radius)/(dx)));
					  (*n1)(dof)=0.5*(1.0-std::tanh((r-radius)/(0.4)));
				  }
			  }
		  }
	  }
  }
}

