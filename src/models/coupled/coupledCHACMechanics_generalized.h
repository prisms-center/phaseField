//Matrix Free implementation of coupled Cahn-Hilliard, Allen-Cahn and Mechanics formulation 
#ifndef CHACMECHANICS_H
#define CHACMECHANICS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

//material models
#include "../mechanics/computeStress.h"

// BC object declaration
template <int dim>
class varBCs
{
	public:
	//varBCs();
	std::vector<std::string> var_BC_type;
	std::vector<double> var_BC_val;
};

template <int dim>
class CoupledCHACMechanicsProblem: public MatrixFreePDE<dim>
{
 public: 
  CoupledCHACMechanicsProblem();

  void shiftConcentration();

  void setBCs();

  void inputBCs(int var, int component, std::string BC_type_dim1_min, double BC_value_dim1_min,
  			std::string BC_type_dim1_max, double BC_value_dim1_max, std::string BC_type_dim2_min, double BC_value_dim2_min,
  			std::string BC_type_dim2_max, double BC_value_dim2_max,std::string BC_type_dim3_min, double BC_value_dim3_min,
  			std::string BC_type_dim3_max, double BC_value_dim3_max);

  void inputBCs(int var, int component, std::string BC_type_dim1_min, double BC_value_dim1_min,
			std::string BC_type_dim1_max, double BC_value_dim1_max, std::string BC_type_dim2_min, double BC_value_dim2_min,
			std::string BC_type_dim2_max, double BC_value_dim2_max);

  void inputBCs(int var, int component, std::string BC_type_dim1_min, double BC_value_dim1_min,
  			std::string BC_type_dim1_max, double BC_value_dim1_max);

  void inputBCs(int var, int component, std::string BC_type, double BC_value);

  // Boundary condition object
    std::vector<varBCs<dim>> BC_list;

 private:


  // Elasticity matrix variables
  Table<2, double> CIJ_table;
  Table<2, double> CIJ_alpha_table;
  Table<2, double> CIJ_beta_table;
  const static unsigned int CIJ_tensor_size = 2*dim-1+dim/3;
  dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > CIJ, CIJ_alpha, CIJ_beta, CIJ_diff;

  bool c_dependent_misfit;

  Threads::Mutex assembler_lock;

  // Variables needed to calculate the LHS
  std::vector<variable_info<dim>> varInfoListRHS;
  std::vector<variable_info<dim>> resInfoListRHS;

  // Variables needed to calculate the LHS
  unsigned int num_var_LHS;
  std::vector<variable_info<dim>> varInfoListLHS;

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

  void markBoundaries(int field_index);

  void getEnergy(const MatrixFree<dim,double> &data,
    				    std::vector<vectorType*> &dst,
    				    const std::vector<vectorType*> &src,
    				    const std::pair<unsigned int,unsigned int> &cell_range);


  void residualRHS(const std::vector<modelVariable<dim>> & modelVarList,
		  	  	  	  	  	  	  	  	  	  	  	  	  std::vector<modelResidual<dim>> & modelResidualsList) const;

  void residualLHS(const std::vector<modelVariable<dim>> & modelVarList,
  		  	  	  	  	  	  	  	  	  	  	  	  	  modelResidual<dim> & modelRes) const;

  void energyDensity(const std::vector<modelVariable<dim>> & modelVarList, const dealii::VectorizedArray<double> & JxW_value);

};

//constructor
template <int dim>
CoupledCHACMechanicsProblem<dim>::CoupledCHACMechanicsProblem(): MatrixFreePDE<dim>(),
  CIJ_table(CIJ_tensor_size,CIJ_tensor_size), CIJ_alpha_table(CIJ_tensor_size,CIJ_tensor_size), CIJ_beta_table(CIJ_tensor_size,CIJ_tensor_size)
{
  //initialize elasticity matrix
#if defined(MaterialModelV) && defined(MaterialConstantsV)
	if (n_dependent_stiffness == true){
		double materialConstants[]=MaterialConstantsV;
		getCIJMatrix<dim>(MaterialModelV, materialConstants, CIJ_alpha_table, this->pcout);

		double materialConstantsBeta[]=MaterialConstantsBetaV;
		getCIJMatrix<dim>(MaterialModelBetaV, materialConstantsBeta, CIJ_beta_table, this->pcout);

		for (unsigned int i=0; i<CIJ_tensor_size; i++){
			for (unsigned int j=0; j<CIJ_tensor_size; j++){
				CIJ_beta[i][j] =  CIJ_beta_table(i,j);
				CIJ_alpha[i][j] =  CIJ_alpha_table(i,j);
				CIJ_diff[i][j] =  CIJ_beta_table(i,j) - CIJ_alpha_table(i,j);
			}
		}
	}
	else{
		double materialConstants[]=MaterialConstantsV;
		getCIJMatrix<dim>(MaterialModelV, materialConstants, CIJ_table, this->pcout);

		for (unsigned int i=0; i<CIJ_tensor_size; i++){
			for (unsigned int j=0; j<CIJ_tensor_size; j++){
				CIJ[i][j] =  CIJ_table(i,j);
			}
		}
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

// Load variable information for calculating the RHS
varInfoListRHS.reserve(num_var);
unsigned int field_number = 0;
unsigned int scalar_var_index = 0;
unsigned int vector_var_index = 0;
for (unsigned int i=0; i<num_var; i++){
	variable_info<dim> varInfo;
	if (need_value[i] or need_gradient[i] or need_hessian[i]){
		varInfo.global_var_index = i;
		varInfo.global_field_index = field_number;
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

	if (var_type[i] == "SCALAR"){
		field_number++;
	}
	else {
		field_number+=dim;
	}
}

variable_info<dim> resInfoLHS;
	for (unsigned int i=0; i<num_var_LHS; i++){
		if (MatrixFreePDE<dim>::currentFieldIndex == varInfoListLHS[i].global_field_index){
			resInfoLHS = varInfoListLHS[i];
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
field_number = 0;
scalar_var_index = 0;
vector_var_index = 0;
for (unsigned int i=0; i<num_var; i++){
	variable_info<dim> varInfo;
	if (need_value_LHS[i] or need_gradient_LHS[i] or need_hessian_LHS[i]){
		varInfo.global_var_index = i;
		varInfo.global_field_index = field_number;
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

	if (var_type[i] == "SCALAR"){
		field_number++;
	}
	else {
		field_number+=dim;
	}
}

}

#endif
