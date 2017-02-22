//constructor and destructor for matrixFreePDE class

#ifndef MATRIXFREEPDE_MATRIXFREE_H
#define MATRIXFREEPDE_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "init.cc"
#include "reinit.cc"
#include "initForTests.cc"
#include "refine.cc"
#include "invM.cc"
#include "computeLHS.cc"
#include "computeRHS.cc"
#include "modifyFields.cc"
#include "solve.cc"
#include "solveIncrement.cc"
#include "outputResults.cc"
#include "markBoundaries.cc"
#include "boundaryConditions.cc"
#include "initialConditions.cc"
#include "utilities.cc"
#include "calcFreeEnergy.cc"
#include "integrate_and_shift_field.cc"
#include "getOutputTimeSteps.cc"
#include "buildFields.cc"

 //constructor
 template <int dim>
 MatrixFreePDE<dim>::MatrixFreePDE ()
 :
 Subscriptor(),
 triangulation (MPI_COMM_WORLD),
 isTimeDependentBVP(false),
 isEllipticBVP(false),
 dtValue(0.0),
 currentTime(0.0),
 finalTime(0.0),
 currentIncrement(0),
 totalIncrements(1),
 pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
 computing_timer (pcout, TimerOutput::summary, TimerOutput::wall_times)
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


	#ifdef need_val_LHS
	{bool temp[] = need_val_LHS;
	vectorLoad(temp, sizeof(temp), need_value_LHS);}
	#else
	for (unsigned int i=0; i<num_var; i++)
		need_value_LHS.push_back(false);
	#endif
	#ifdef need_grad_LHS
	{bool temp[] = need_grad_LHS;
	vectorLoad(temp, sizeof(temp), need_gradient_LHS);}
	#else
	for (unsigned int i=0; i<num_var; i++)
		need_gradient_LHS.push_back(false);
	#endif
	#ifdef need_hess_LHS
	{bool temp[] = need_hess_LHS;
	vectorLoad(temp, sizeof(temp), need_hessian_LHS);}
	#else
	for (unsigned int i=0; i<num_var; i++)
		need_hessian_LHS.push_back(false);
	#endif
	#ifdef need_val_residual_LHS
	{bool temp[] = need_val_residual_LHS;
	vectorLoad(temp, sizeof(temp), value_residual_LHS);}
	#else
	for (unsigned int i=0; i<num_var; i++)
		value_residual_LHS.push_back(false);
	#endif
	#ifdef need_grad_residual_LHS
	{bool temp[] = need_grad_residual_LHS;
	vectorLoad(temp, sizeof(temp), gradient_residual_LHS);}
	#else
	for (unsigned int i=0; i<num_var; i++)
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

 //destructor
 template <int dim>
 MatrixFreePDE<dim>::~MatrixFreePDE ()
 {
   matrixFreeObject.clear();
   for(unsigned int iter=0; iter<fields.size(); iter++){
     delete soltransSet[iter];
     delete locally_relevant_dofsSet[iter];
     delete constraintsDirichletSet[iter];
     delete dofHandlersSet[iter];
     delete FESet[iter];
     delete solutionSet[iter];
     delete residualSet[iter];
   } 
 }

#endif
