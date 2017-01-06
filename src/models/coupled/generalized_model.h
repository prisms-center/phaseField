//Matrix Free implementation of general system of PDEs

//this source file is temporarily treated as a header file until library packaging scheme is finalized

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
class generalizedProblem: public MatrixFreePDE<dim>
{
 public: 
  generalizedProblem();

  void shiftConcentration();

  void setBCs();

  void buildFields();

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
    std::vector<varBCs<dim> > BC_list;

 private:

    // Declare vectors to put the information about the variables into
    std::vector<std::string> var_name;
    std::vector<std::string> var_type;
    std::vector<std::string> var_eq_type;

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

  // Elasticity matrix variables
  const static unsigned int CIJ_tensor_size = 2*dim-1+dim/3;
  std::vector<dealii::Tensor<2, CIJ_tensor_size, dealii::VectorizedArray<double> > > CIJ_list;

  bool c_dependent_misfit;

  Threads::Mutex assembler_lock;

  // Variables needed to calculate the LHS
  std::vector<variable_info<dim> > varInfoListRHS;
  std::vector<variable_info<dim> > resInfoListRHS;

  // Variables needed to calculate the LHS
  unsigned int num_var_LHS;
  std::vector<variable_info<dim> > varInfoListLHS;

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

  void getEnergy(const MatrixFree<dim,double> &data,
    				    std::vector<vectorType*> &dst,
    				    const std::vector<vectorType*> &src,
    				    const std::pair<unsigned int,unsigned int> &cell_range);


  void residualRHS(const std::vector<modelVariable<dim> > & modelVarList,
		  	  	  	  	  	  	  	  	  	  	  	  	  std::vector<modelResidual<dim> > & modelResidualsList,
														  dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

  void residualLHS(const std::vector<modelVariable<dim> > & modelVarList,
  		  	  	  	  	  	  	  	  	  	  	  	  	  modelResidual<dim> & modelRes,
														  dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const;

  void energyDensity(const std::vector<modelVariable<dim> > & modelVarList, const dealii::VectorizedArray<double> & JxW_value,
		  	  	  	  	  	  	  	  	  	  	  	  	  dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc);

  //AMR methods
  void adaptiveRefine(unsigned int currentIncrement);
  void adaptiveRefineCriterion();

};

