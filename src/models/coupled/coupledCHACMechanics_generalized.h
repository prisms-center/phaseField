//Matrix Free implementation of coupled Cahn-Hilliard, Allen-Cahn and Mechanics formulation 
#ifndef CHACMECHANICS_H
#define CHACMECHANICS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

//material models
#include "../mechanics/computeStress.h"

template <int dim>
class CoupledCHACMechanicsProblem: public MatrixFreePDE<dim>
{
 public: 
  CoupledCHACMechanicsProblem();

  void shiftConcentration();

 private:
  //elasticity matrix
  Table<2, double> CIJ;
  Table<2, double> CIJ_alpha;
  Table<2, double> CIJ_beta;
  Table<2, double> CIJ_diff;

  dealii::Tensor<2, 2*dim-1+dim/3, dealii::VectorizedArray<double> > CIJ_alpha_tensor, CIJ_beta_tensor;

  //RHS implementation for explicit solve
  void getRHS(const MatrixFree<dim,double> &data, 
	      std::vector<vectorType*> &dst, 
	      const std::vector<vectorType*> &src,
	      const std::pair<unsigned int,unsigned int> &cell_range) const;
    
  //LHS implementation for implicit solve 
  void  getLHS(const MatrixFree<dim,double> &data, 
	       vectorType &dst, 
	       const vectorType &src,
	       const std::pair<unsigned int,unsigned int> &cell_range) const;

  //method to apply initial conditions
  void applyInitialConditions();
 
  //methods to apply dirichlet BC's on displacement
  void applyDirichletBCs();

  // method to modify the fields for nucleation
  void modifySolutionFields();

  void computeIntegral(double& integratedField);

  void markBoundaries();

  bool c_dependent_misfit;

  void getEnergy(const MatrixFree<dim,double> &data,
    				    std::vector<vectorType*> &dst,
    				    const std::vector<vectorType*> &src,
    				    const std::pair<unsigned int,unsigned int> &cell_range);
  Threads::Mutex assembler_lock;

  void residualRHS(const std::vector<modelVariable<dim>> & modelVariablesList,
		  	  	  	  	  	  	  	  	  	  	  	  	  std::vector<modelResidual<dim>> & modelResidualsList) const;

  void residualLHS(const std::vector<modelVariable<dim>> & modelVariablesList,
  		  	  	  	  	  	  	  	  	  	  	  	  	  std::vector<modelResidual<dim>> & modelResidualsList) const;

  unsigned int num_var_LHS;
  unsigned int num_res_LHS;
  std::vector<unsigned int> LHS_indices;
  std::vector<unsigned int> LHS_res_indices;
  std::vector<unsigned int> LHS_var_indices_of_res;

};

//constructor
template <int dim>
CoupledCHACMechanicsProblem<dim>::CoupledCHACMechanicsProblem(): MatrixFreePDE<dim>(),
  CIJ(2*dim-1+dim/3,2*dim-1+dim/3), CIJ_alpha(2*dim-1+dim/3,2*dim-1+dim/3), CIJ_beta(2*dim-1+dim/3,2*dim-1+dim/3), CIJ_diff(2*dim-1+dim/3,2*dim-1+dim/3)
{
  //initialize elasticity matrix
#if defined(MaterialModelV) && defined(MaterialConstantsV)
	if (n_dependent_stiffness == true){
		double materialConstants[]=MaterialConstantsV;
		getCIJMatrix<dim>(MaterialModelV, materialConstants, CIJ_alpha, this->pcout);

		double materialConstantsBeta[]=MaterialConstantsBetaV;
		getCIJMatrix<dim>(MaterialModelBetaV, materialConstantsBeta, CIJ_beta, this->pcout);

		for (unsigned int i=0; i<2*dim-1+dim/3; i++){
			for (unsigned int j=0; j<2*dim-1+dim/3; j++){
				CIJ_beta_tensor[i][j] =  CIJ_beta(i,j);
				CIJ_alpha_tensor[i][j] =  CIJ_alpha(i,j);
				CIJ_diff(i,j) = CIJ_beta(i,j) - CIJ_alpha(i,j);
			}
		}
	}
	else{
		double materialConstants[]=MaterialConstantsV;
		getCIJMatrix<dim>(MaterialModelV, materialConstants, CIJ, this->pcout);
	}

#else
#error Compile ERROR: missing material property variable: MaterialModelV, MaterialConstantsV
#endif

c_dependent_misfit = false;
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		if ((std::abs(sfts_linear1[i][j])>1.0e-12)||(std::abs(sfts_linear2[i][j])>1.0e-12)||(std::abs(sfts_linear3[i][j])>1.0e-12)){
			c_dependent_misfit = true;
		}
	}
}

// If interpolation functions for the strain aren't specifically defined, use the general interpolation functions
#ifndef h1strainV
	#define h1strainV h1V
#endif
#ifndef h2strainV
	#define h2strainV h2V
#endif
#ifndef h3strainV
	#define h3strainV h3V
#endif

#ifndef hn1strainV
	#define hn1strainV hn1V
#endif
#ifndef hn2strainV
	#define hn2strainV hn2V
#endif
#ifndef hn3strainV
	#define hn3strainV hn3V
#endif

// If the Landau energy terms aren't defined, set them to zero
#ifndef W
	#define W 0.0
#endif
#ifndef fbarrierV
	#define fbarrierV 0.0
#endif

// If nucleation isn't specifically turned on, set nucleation_occurs to false
#ifndef nucleation_occurs
	#define nucleation_occurs false
#endif

num_var_LHS = 0;
num_res_LHS = 0;
for (unsigned int i=0; i<num_var; i++){
	if (need_value_LHS[i] or need_gradient_LHS[i] or need_hessian_LHS[i]){
		num_var_LHS++;
		LHS_indices.push_back(i);
	}
	if (value_residual_LHS[i] or gradient_residual_LHS[i]){
		num_res_LHS++;
		LHS_res_indices.push_back(i);
		LHS_var_indices_of_res.push_back(num_var_LHS-1);
	}
}

}

template <int dim>
void CoupledCHACMechanicsProblem<dim>::getRHS(const MatrixFree<dim,double> &data,
					       std::vector<vectorType*> &dst, 
					       const std::vector<vectorType*> &src,
					       const std::pair<unsigned int,unsigned int> &cell_range) const{


  //initialize FEEvaulation objects
  std::vector<typeScalar> scalar_vars;
  std::vector<typeVector> vector_vars;
  std::vector<bool> is_scalar_var;
  std::vector<unsigned int> scalar_or_vector_index;
  std::vector<unsigned int> field_index;
  unsigned int field_number = 0;

  for (unsigned int i=0; i<num_var; i++){
	  if (var_type[i] == "SCALAR"){
		  typeScalar var(data, i);
		  scalar_vars.push_back(var);
		  is_scalar_var.push_back(true);
		  scalar_or_vector_index.push_back(scalar_vars.size()-1);
		  field_index.push_back(field_number);
		  field_number++;
	  }
	  else {
		  typeVector var(data, i);
		  vector_vars.push_back(var);
		  is_scalar_var.push_back(false);
		  scalar_or_vector_index.push_back(vector_vars.size()-1);
		  field_index.push_back(field_number);
		  field_number+=dim;
	  }
  }

  std::vector<modelVariable<dim> > modelVariablesList;
  std::vector<modelResidual<dim> > modelResidualsList;
  modelVariablesList.reserve(num_var);
  modelResidualsList.reserve(num_var);

  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){

	  // Initialize, read DOFs, and set evaulation flags for each variable
	  for (unsigned int i=0; i<num_var; i++){
		  if (is_scalar_var[i] == true) {
			  scalar_vars[scalar_or_vector_index[i]].reinit(cell);
			  scalar_vars[scalar_or_vector_index[i]].read_dof_values_plain(*src[field_index[i]]);
			  scalar_vars[scalar_or_vector_index[i]].evaluate(need_value[i], need_gradient[i], need_hessian[i]);
		  }
		  else {
			  vector_vars[scalar_or_vector_index[i]].reinit(cell);
			  vector_vars[scalar_or_vector_index[i]].read_dof_values_plain(*src[field_index[i]]);
			  vector_vars[scalar_or_vector_index[i]].evaluate(need_value[i], need_gradient[i], need_hessian[i]);
		  }
	  }

	  unsigned int n_q_points = scalar_vars[0].n_q_points;

	  //loop over quadrature points
	  for (unsigned int q=0; q<n_q_points; ++q){

		  for (unsigned int i=0; i<num_var; i++){
			  if (is_scalar_var[i] == true) {
				  if (need_value[i] == true){
					  modelVariablesList[i].scalarValue = scalar_vars[scalar_or_vector_index[i]].get_value(q);
				  }
				  if (need_gradient[i] == true){
					  modelVariablesList[i].scalarGrad = scalar_vars[scalar_or_vector_index[i]].get_gradient(q);
				  }
				  if (need_hessian[i] == true){
					  modelVariablesList[i].scalarHess = scalar_vars[scalar_or_vector_index[i]].get_hessian(q);
				  }
			  }
			  else {
				  if (need_value[i] == true){
					  modelVariablesList[i].vectorValue = vector_vars[scalar_or_vector_index[i]].get_value(q);
				  }
				  if (need_gradient[i] == true){
					  modelVariablesList[i].vectorGrad = vector_vars[scalar_or_vector_index[i]].get_gradient(q);
				  }
				  if (need_hessian[i] == true){
					  modelVariablesList[i].vectorHess = vector_vars[scalar_or_vector_index[i]].get_hessian(q);
				  }
			  }
		  }

		  // Calculate the residuals
		  residualRHS(modelVariablesList,modelResidualsList);
  
		  // Submit values
		  for (unsigned int i=0; i<num_var; i++){
			  if (is_scalar_var[i] == true) {
				  if (value_residual[i] == true){
					  scalar_vars[scalar_or_vector_index[i]].submit_value(modelResidualsList[i].scalarValueResidual,q);
				  }
      			  if (gradient_residual[i] == true){
      				  scalar_vars[scalar_or_vector_index[i]].submit_gradient(modelResidualsList[i].scalarGradResidual,q);
      			  }
      		  }
      		  else {
      			  if (value_residual[i] == true){
      				  vector_vars[scalar_or_vector_index[i]].submit_value(modelResidualsList[i].vectorValueResidual,q);
      			  }
      			  if (gradient_residual[i] == true){
      				  vector_vars[scalar_or_vector_index[i]].submit_gradient(modelResidualsList[i].vectorGradResidual,q);
      			  }
      		  }
      	  }

	  }

	  for (unsigned int i=0; i<num_var; i++){
		  if (is_scalar_var[i] == true) {
			  scalar_vars[scalar_or_vector_index[i]].integrate(value_residual[i], gradient_residual[i]);
			  scalar_vars[scalar_or_vector_index[i]].distribute_local_to_global(*dst[field_index[i]]);
		  }
		  else {
			  vector_vars[scalar_or_vector_index[i]].integrate(value_residual[i], gradient_residual[i]);
			  vector_vars[scalar_or_vector_index[i]].distribute_local_to_global(*dst[field_index[i]]);
		  }
	  }
  }
}

template <int dim>
void  CoupledCHACMechanicsProblem<dim>::getLHS(const MatrixFree<dim,double> &data, 
					       vectorType &dst, 
					       const vectorType &src,
					       const std::pair<unsigned int,unsigned int> &cell_range) const{

	//initialize FEEvaulation objects
	std::vector<typeScalar> scalar_vars;
	std::vector<typeVector> vector_vars;
	std::vector<bool> is_scalar_var;
	std::vector<unsigned int> scalar_or_vector_index;
	std::vector<unsigned int> field_index;
	unsigned int field_number = 0;

//	  typeScalar var1(data,1);
//	  scalar_vars.push_back(var1);
//	  is_scalar_var.push_back(true);
//	  scalar_or_vector_index.push_back(0);
//	  field_index.push_back(1);
//	  typeScalar var2(data,2);
//	  scalar_vars.push_back(var2);
//	  is_scalar_var.push_back(true);
//	  scalar_or_vector_index.push_back(1);
//	  field_index.push_back(2);
//	  typeScalar var3(data,3);
//	  scalar_vars.push_back(var3);
//	  is_scalar_var.push_back(true);
//	  scalar_or_vector_index.push_back(2);
//	  field_index.push_back(4);
//	  typeVector var4(data,4);
//	  vector_vars.push_back(var4);
//	  is_scalar_var.push_back(true);
//	  scalar_or_vector_index.push_back(0);
//	  field_index.push_back(5);

	  for (unsigned int i=0; i<num_var_LHS; i++){
	  		if (var_type[LHS_indices[i]] == "SCALAR"){
	  			typeScalar var(data, LHS_indices[i]);
	  			scalar_vars.push_back(var);
	  			is_scalar_var.push_back(true);
	  			scalar_or_vector_index.push_back(scalar_vars.size()-1);
	  			field_index.push_back(field_number);
	  			field_number++;
	  		}
	  		  else {
	  			  typeVector var(data, LHS_indices[i]);
	  			  vector_vars.push_back(var);
	  			  is_scalar_var.push_back(false);
	  			  scalar_or_vector_index.push_back(vector_vars.size()-1);
	  			  field_index.push_back(field_number);
	  			  field_number+=dim;
	  		  }
	  	  }

	  std::vector<modelVariable<dim> > modelVariablesList;
	  std::vector<modelResidual<dim> > modelResidualsList;
	  modelVariablesList.reserve(num_var_LHS);
	  modelResidualsList.reserve(num_res_LHS);

	  //loop over cells

		    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){

		  	  // Initialize, read DOFs, and set evaulation flags for each variable
		  	  for (unsigned int i=0; i<num_var_LHS; i++){
		  		  if (is_scalar_var[i] == true) {
		  			  scalar_vars[scalar_or_vector_index[i]].reinit(cell);
		  			  if ( (value_residual_LHS[LHS_indices[i]] == true) or (gradient_residual_LHS[LHS_indices[i]] == true) ){
		  				  scalar_vars[scalar_or_vector_index[i]].read_dof_values_plain(src);  // need to generalize, I think
		  			  }
		  			  else{
		  				  scalar_vars[scalar_or_vector_index[i]].read_dof_values_plain(*MatrixFreePDE<dim>::solutionSet[LHS_indices[i]]);
		  			  }
		  			  scalar_vars[scalar_or_vector_index[i]].evaluate(need_value_LHS[LHS_indices[i]], need_gradient_LHS[LHS_indices[i]], need_hessian_LHS[LHS_indices[i]]);
		  		  }
		  		  else {
		  			  vector_vars[scalar_or_vector_index[i]].reinit(cell);
		  			  if ( (value_residual_LHS[LHS_indices[i]] == true) or (gradient_residual_LHS[LHS_indices[i]] == true) ){
		  				  vector_vars[scalar_or_vector_index[i]].read_dof_values_plain(src);
		  			  }
		  			  else {
		  				  vector_vars[scalar_or_vector_index[i]].read_dof_values_plain(*MatrixFreePDE<dim>::solutionSet[field_index[i]]);
		  			  }
		  			  vector_vars[scalar_or_vector_index[i]].evaluate(need_value_LHS[LHS_indices[i]], need_gradient_LHS[LHS_indices[i]], need_hessian_LHS[LHS_indices[i]]);
		  		  }
		  	  }

	    //loop over quadrature points
	    for (unsigned int q=0; q<scalar_vars[0].n_q_points; ++q){
//	      //u
//	      vectorgradType ux = vector_vars[0].get_gradient(q);
//	      vectorgradType Rux;
//
//	      // n
//	      scalarvalueType n1, n2, n3;
//	      n1 = scalar_vars[0].get_value(q);
//		  #if num_sop>1
//			  n2 = scalar_vars[1].get_value(q);
//		  #else
//			  n2 = constV(0.0);
//		  #endif
//		  #if num_sop>2
//			  n3 = scalar_vars[2].get_value(q);
//		  #else
//			  n3 = constV(0.0);
//		  #endif

			  for (unsigned int i=0; i<num_var_LHS; i++){
				  if (is_scalar_var[i] == true) {
					  if (need_value_LHS[LHS_indices[i]] == true){
						  modelVariablesList[i].scalarValue = scalar_vars[scalar_or_vector_index[i]].get_value(q);
					  }
					  if (need_gradient_LHS[LHS_indices[i]] == true){
						  modelVariablesList[i].scalarGrad = scalar_vars[scalar_or_vector_index[i]].get_gradient(q);
					  }
					  if (need_hessian_LHS[LHS_indices[i]] == true){
						  modelVariablesList[i].scalarHess = scalar_vars[scalar_or_vector_index[i]].get_hessian(q);
					  }
				  }
				  else {
					  if (need_value_LHS[LHS_indices[i]] == true){
						  modelVariablesList[i].vectorValue = vector_vars[scalar_or_vector_index[i]].get_value(q);
					  }
					  if (need_gradient_LHS[LHS_indices[i]] == true){
						  modelVariablesList[i].vectorGrad = vector_vars[scalar_or_vector_index[i]].get_gradient(q);
					  }
					  if (need_hessian_LHS[LHS_indices[i]] == true){
						  modelVariablesList[i].vectorHess = vector_vars[scalar_or_vector_index[i]].get_hessian(q);
					  }
				  }
			  }

			  //n1
			  scalarvalueType n1 = modelVariablesList[0].scalarValue;

			  //n2
			  scalarvalueType n2 = modelVariablesList[1].scalarValue;


			  //n3
			  scalarvalueType n3 = modelVariablesList[2].scalarValue;

			  //u
			  vectorgradType ux = modelVariablesList[3].vectorGrad;
			  vectorgradType Rux;

			  // Calculate the residuals
			  residualLHS(modelVariablesList,modelResidualsList);

			  // Submit values
			  for (unsigned int i=0; i<num_res_LHS; i++){
				  if (var_type[LHS_res_indices[i]] == "SCALAR"){
					  if (value_residual[LHS_res_indices[i]] == true){
						  scalar_vars[scalar_or_vector_index[LHS_var_indices_of_res[i]]].submit_value(modelResidualsList[scalar_or_vector_index[i]].scalarValueResidual,q);
					  }
					  if (gradient_residual[LHS_res_indices[i]] == true){
						  scalar_vars[scalar_or_vector_index[LHS_var_indices_of_res[i]]].submit_gradient(modelResidualsList[scalar_or_vector_index[i]].scalarGradResidual,q);
					  }
				  }
				  else {
					  if (value_residual[LHS_res_indices[i]] == true){
						  vector_vars[scalar_or_vector_index[LHS_var_indices_of_res[i]]].submit_value(modelResidualsList[scalar_or_vector_index[i]].vectorValueResidual,q);
					  }
					  if (gradient_residual[LHS_res_indices[i]] == true){
						  vector_vars[scalar_or_vector_index[LHS_var_indices_of_res[i]]].submit_gradient(modelResidualsList[scalar_or_vector_index[i]].vectorGradResidual,q);
					  }
				  }
			  }
	    }

	    //integrate
	  for (unsigned int i=0; i<num_res_LHS; i++){
		  if (var_type[LHS_res_indices[i]] == "SCALAR") {
			  scalar_vars[scalar_or_vector_index[LHS_var_indices_of_res[i]]].integrate(value_residual[LHS_res_indices[i]], gradient_residual[LHS_res_indices[i]]);
			  scalar_vars[scalar_or_vector_index[LHS_var_indices_of_res[i]]].distribute_local_to_global(dst);
		  }
		  else {
			  vector_vars[scalar_or_vector_index[LHS_var_indices_of_res[i]]].integrate(value_residual[LHS_res_indices[i]], gradient_residual[LHS_res_indices[i]]);
			  vector_vars[scalar_or_vector_index[LHS_var_indices_of_res[i]]].distribute_local_to_global(dst);
		  }
	  }
	  }

//	//initialize FEEvaulation objects
//	std::vector<typeScalar> scalar_vars;
//	std::vector<typeVector> vector_vars;
//	std::vector<bool> is_scalar_var;
//	std::vector<unsigned int> scalar_or_vector_index;
//	std::vector<unsigned int> field_index;
//	unsigned int field_number = 0;
//
//	for (unsigned int i=0; i<num_var_LHS; i++){
//		if (var_type[LHS_indices[i]] == "SCALAR"){
//			typeScalar var(data, LHS_indices[i]);
//			scalar_vars.push_back(var);
//			is_scalar_var.push_back(true);
//			scalar_or_vector_index.push_back(scalar_vars.size()-1);
//			field_index.push_back(field_number);
//			field_number++;
//		}
//		  else {
//			  typeVector var(data, LHS_indices[i]);
//			  vector_vars.push_back(var);
//			  is_scalar_var.push_back(false);
//			  scalar_or_vector_index.push_back(vector_vars.size()-1);
//			  field_index.push_back(field_number);
//			  field_number+=dim;
//		  }
//	  }
//
//	  std::vector<modelVariable<dim> > modelVariablesList;
//	  std::vector<modelResidual<dim> > modelResidualsList;
//	  modelVariablesList.reserve(num_var_LHS);
//	  modelResidualsList.reserve(num_res_LHS);
//
//	  typeScalar var(data,1);
//	  scalar_Vars.push_back(var);
//	  is_scalar_var.push_back(true);
//	  scalar_or_vector_index.push_back(0);
//	  field_index.push_back(1);
//	  typeScalar var(data,2);
//	  scalar_Vars.push_back(var);
//	  is_scalar_var.push_back(true);
//	  scalar_or_vector_index.push_back(1);
//	  field_index.push_back(2);
//	  typeScalar var(data,3);
//	  scalar_Vars.push_back(var);
//	  is_scalar_var.push_back(true);
//	  scalar_or_vector_index.push_back(2);
//	  field_index.push_back(4);
//	  typeVector var(data,4);
//	  vector_Vars.push_back(var);
//	  is_scalar_var.push_back(true);
//	  scalar_or_vector_index.push_back(0);
//	  field_index.push_back(5);
//
//  //loop over cells
//  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
//
//	  // Initialize, read DOFs, and set evaulation flags for each variable
//	  for (unsigned int i=0; i<num_var_LHS; i++){
//		  if (is_scalar_var[i] == true) {
//			  scalar_vars[scalar_or_vector_index[i]].reinit(cell);
//			  if ( (value_residual_LHS[LHS_indices[i]] == true) or (gradient_residual_LHS[LHS_indices[i]] == true) ){
//				  scalar_vars[scalar_or_vector_index[i]].read_dof_values_plain(src);  // need to generalize, I think
//			  }
//			  else{
//				  scalar_vars[scalar_or_vector_index[i]].read_dof_values_plain(*MatrixFreePDE<dim>::solutionSet[LHS_indices[i]]);
//			  }
//			  scalar_vars[scalar_or_vector_index[i]].evaluate(need_value_LHS[LHS_indices[i]], need_gradient_LHS[LHS_indices[i]], need_hessian_LHS[LHS_indices[i]]);
//		  }
//		  else {
//			  vector_vars[scalar_or_vector_index[i]].reinit(cell);
//			  if ( (value_residual_LHS[LHS_indices[i]] == true) or (gradient_residual_LHS[LHS_indices[i]] == true) ){
//				  vector_vars[scalar_or_vector_index[i]].read_dof_values_plain(src);
//			  }
//			  else {
//				  vector_vars[scalar_or_vector_index[i]].read_dof_values_plain(*MatrixFreePDE<dim>::solutionSet[field_index[i]]);
//			  }
//			  vector_vars[scalar_or_vector_index[i]].evaluate(need_value_LHS[LHS_indices[i]], need_gradient_LHS[LHS_indices[i]], need_hessian_LHS[LHS_indices[i]]);
//		  }
//	  }
//
//	  unsigned int n_q_points = scalar_vars[0].n_q_points;
//
//
//	  //loop over quadrature points
//	  for (unsigned int q=0; q<n_q_points; ++q){
//		  for (unsigned int i=0; i<num_var_LHS; i++){
//			  if (is_scalar_var[i] == true) {
//				  if (need_value_LHS[LHS_indices[i]] == true){
//					  modelVariablesList[i].scalarValue = scalar_vars[scalar_or_vector_index[i]].get_value(q);
//				  }
//				  if (need_gradient_LHS[LHS_indices[i]] == true){
//					  modelVariablesList[i].scalarGrad = scalar_vars[scalar_or_vector_index[i]].get_gradient(q);
//				  }
//				  if (need_hessian_LHS[LHS_indices[i]] == true){
//					  modelVariablesList[i].scalarHess = scalar_vars[scalar_or_vector_index[i]].get_hessian(q);
//				  }
//			  }
//			  else {
//				  if (need_value_LHS[LHS_indices[i]] == true){
//					  modelVariablesList[i].vectorValue = vector_vars[scalar_or_vector_index[i]].get_value(q);
//				  }
//				  if (need_gradient_LHS[LHS_indices[i]] == true){
//					  modelVariablesList[i].vectorGrad = vector_vars[scalar_or_vector_index[i]].get_gradient(q);
//				  }
//				  if (need_hessian_LHS[LHS_indices[i]] == true){
//					  modelVariablesList[i].vectorHess = vector_vars[scalar_or_vector_index[i]].get_hessian(q);
//				  }
//			  }
//		  }
//
//		  // Calculate the residuals
//		  residualLHS(modelVariablesList,modelResidualsList);
//
//		  // Submit values
//		  for (unsigned int i=0; i<num_res_LHS; i++){
//			  if (var_type[LHS_res_indices[i]] == "SCALAR"){
//				  if (value_residual[LHS_res_indices[i]] == true){
//					  scalar_vars[LHS_var_indices_of_res[i]].submit_value(modelResidualsList[LHS_res_indices[i]].scalarValueResidual,q);
//				  }
//				  if (gradient_residual[LHS_res_indices[i]] == true){
//					  scalar_vars[LHS_var_indices_of_res[i]].submit_gradient(modelResidualsList[LHS_res_indices[i]].scalarGradResidual,q);
//				  }
//			  }
//			  else {
//				  if (value_residual[LHS_res_indices[i]] == true){
//					  vector_vars[LHS_var_indices_of_res[i]].submit_value(modelResidualsList[LHS_res_indices[i]].vectorValueResidual,q);
//				  }
//				  if (gradient_residual[LHS_res_indices[i]] == true){
//					  vector_vars[LHS_var_indices_of_res[i]].submit_gradient(modelResidualsList[LHS_res_indices[i]].vectorGradResidual,q);
//				  }
//			  }
//		  }
//
//	  }
//
//	  for (unsigned int i=0; i<num_res_LHS; i++){
//		  if (var_type[LHS_res_indices[i]] == "SCALAR") {
//			  scalar_vars[LHS_var_indices_of_res[i]].integrate(value_residual[LHS_res_indices[i]], gradient_residual[LHS_res_indices[i]]);
//			  scalar_vars[LHS_var_indices_of_res[i]].distribute_local_to_global(dst);
//		  }
//		  else {
//			  std::cout<<LHS_var_indices_of_res[i] << std::endl;
//			  vector_vars[LHS_var_indices_of_res[i]].integrate(value_residual[LHS_res_indices[i]], gradient_residual[LHS_res_indices[i]]);
//			  vector_vars[LHS_var_indices_of_res[i]].distribute_local_to_global(dst);
//		  }
//	  }
//
//  }
}

// Calculate the free energy
template <int dim>
void  CoupledCHACMechanicsProblem<dim>::getEnergy(const MatrixFree<dim,double> &data,
				    std::vector<vectorType*> &dst,
				    const std::vector<vectorType*> &src,
				    const std::pair<unsigned int,unsigned int> &cell_range) {

	//initialize fields
	  typeScalar cVals(data, 0);

	  typeScalar n1Vals(data,1);
	  #if num_sop>2
	  typeScalar n3Vals(data,3);
	  #endif
	  #if num_sop>1
	  typeScalar n2Vals(data,2);
	  #endif

	  typeVector uVals(data, num_sop+1);

	  //loop over cells
	  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
	    //initialize c field
	    cVals.reinit(cell); cVals.read_dof_values_plain(*src[0]); cVals.evaluate(true, true, false);

	    //initialize n fields
	    n1Vals.reinit(cell); n1Vals.read_dof_values_plain(*src[1]); n1Vals.evaluate(true, true, false);
		#if num_sop>1
	    n2Vals.reinit(cell); n2Vals.read_dof_values_plain(*src[2]); n2Vals.evaluate(true, true, false);
		#endif
		#if num_sop>2
	    n3Vals.reinit(cell); n3Vals.read_dof_values_plain(*src[3]); n3Vals.evaluate(true, true, false);
		#endif

	    //initialize u field
	    uVals.reinit(cell); uVals.read_dof_values_plain(*src[num_sop+1]);
	    uVals.evaluate(false, true, false);

	    dealii::AlignedVector<dealii::VectorizedArray<double> > JxW(cVals.n_q_points);
	    cVals.fill_JxW_values(JxW);

	    //loop over quadrature points
	    for (unsigned int q=0; q<cVals.n_q_points; ++q){
	      //c
	      scalarvalueType c = cVals.get_value(q);
	      scalargradType cx = cVals.get_gradient(q);

	      //n1
	      scalarvalueType n1 = n1Vals.get_value(q);
	      scalargradType n1x = n1Vals.get_gradient(q);

	      //n2
	      scalarvalueType n2, n3;
	      scalargradType n2x, n3x;
		  #if num_sop>1
	    	  n2 = n2Vals.get_value(q);
	    	  n2x = n2Vals.get_gradient(q);
		  #else
	    	  n2 = constV(0.0);
	    	  n2x = constV(0.0)*n1x;
		  #endif

	      //n3
		  #if num_sop>2
	    	  n3 = n3Vals.get_value(q);
	    	  n3x = n3Vals.get_gradient(q);
		  #else
	    	  n3 = constV(0.0);
	    	  n3x = constV(0.0)*n1x;
		  #endif
	      //u
	      vectorgradType ux = uVals.get_gradient(q);

	      scalarvalueType total_energy_density = constV(0.0);

	      scalarvalueType f_chem = (constV(1.0)-(h1V+h2V+h3V))*faV + (h1V+h2V+h3V)*fbV + W*fbarrierV;

	      scalarvalueType f_grad = constV(0.0);

	      for (int i=0; i<dim; i++){
	    	  for (int j=0; j<dim; j++){
	    		  f_grad += constV(0.5*Kn1[i][j])*n1x[i]*n1x[j];
	    	  }
	      }
#if num_sop>1
	      for (int i=0; i<dim; i++){
	    	  for (int j=0; j<dim; j++){
	    		  f_grad += constV(0.5*Kn2[i][j])*n2x[i]*n2x[j];
	    	  }
	      }
#endif
#if num_sop>2
	      for (int i=0; i<dim; i++){
	    	  for (int j=0; j<dim; j++){
	    		  f_grad += constV(0.5*Kn3[i][j])*n3x[i]*n3x[j];
	    	  }
	      }
#endif


	      // Calculate the stress-free transformation strain and its derivatives at the quadrature point
	      dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc, sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

	      for (unsigned int i=0; i<dim; i++){
	    	  for (unsigned int j=0; j<dim; j++){
	    		  // Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c + b_p
	    		  sfts1[i][j] = constV(sfts_linear1[i][j])*c + constV(sfts_const1[i][j]);
	    		  sfts1c[i][j] = constV(sfts_linear1[i][j]);
	    		  sfts1cc[i][j] = constV(0.0);

	    		  // Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c + b_p
	    		  sfts2[i][j] = constV(sfts_linear2[i][j])*c + constV(sfts_const2[i][j]);
	    		  sfts2c[i][j] = constV(sfts_linear1[i][j]);
	    		  sfts2cc[i][j] = constV(0.0);

	    		  // Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c + b_p
	    		  sfts3[i][j] = constV(sfts_linear3[i][j])*c + constV(sfts_const3[i][j]);
	    		  sfts3c[i][j] = constV(sfts_linear3[i][j]);
	    		  sfts3cc[i][j] = constV(0.0);
	    	  }
	      }

	      //compute E2=(E-E0)
	      dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

	      for (unsigned int i=0; i<dim; i++){
	    	  for (unsigned int j=0; j<dim; j++){
	    		  //E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V + sfts2[i][j]*h2V + sfts3[i][j]*h3V);
	    		  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1strainV + sfts2[i][j]*h2strainV + sfts3[i][j]*h3strainV);

	    	  }
	      }

	      //compute stress
	      //S=C*(E-E0)
	      dealii::VectorizedArray<double> CIJ_combined[2*dim-1+dim/3][2*dim-1+dim/3];

	      if (n_dependent_stiffness == true){
	    	  dealii::VectorizedArray<double> sum_hV;
	    	  sum_hV = h1V+h2V+h3V;
	    	  for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	    		  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
	    			  CIJ_combined[i][j] = constV(CIJ_alpha(i,j))*(constV(1.0)-sum_hV) + constV(CIJ_beta(i,j))*sum_hV;
	    		  }
	    	  }
	    	  computeStress<dim>(CIJ_combined, E2, S);
	      }
	      else{
	    	  computeStress<dim>(CIJ, E2, S);
	      }

	      scalarvalueType f_el = constV(0.0);

	      for (unsigned int i=0; i<dim; i++){
	    	  for (unsigned int j=0; j<dim; j++){
	    		  f_el += constV(0.5) * S[i][j]*E2[i][j];
	    	  }
	      }

	      total_energy_density = f_chem + f_grad + f_el;

	      assembler_lock.acquire ();
	      for (unsigned i=0; i<c.n_array_elements;i++){
	    	  // For some reason, some of the values in this loop
	    	  if (c[i] > 1.0e-10){
	    		  this->energy+=total_energy_density[i]*JxW[q][i];
	    		  this->energy_components[0]+= f_chem[i]*JxW[q][i];
	    		  this->energy_components[1]+= f_grad[i]*JxW[q][i];
	    		  this->energy_components[2]+= f_el[i]*JxW[q][i];
	    	  }
	      }
	      assembler_lock.release ();
	    }
	}
}


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
void CoupledCHACMechanicsProblem<dim>::modifySolutionFields()
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

//compute the integral of one of the fields
template <int dim>
void CoupledCHACMechanicsProblem<dim>::computeIntegral(double& integratedField){
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

#endif
