// Class to load in the user input from parameters.h and the variable definition
// part of equations.h

#ifndef INCLUDE_USERINPUTPARAMETERS_H_
#define INCLUDE_USERINPUTPARAMETERS_H_

#include "dealIIheaders.h"
#include "list_of_CIJ.h"
#include "../src/userInputParameters/getCIJMatrix.h"
#include "model_variables.h"
#include "varBCs.h"

template <int dim>
class userInputParameters
{
private:
	// Method to create the list of time steps where the results should be output (called from loadInputParameters)
	std::vector<unsigned int> setOutputTimeSteps(const std::string outputSpacingType, unsigned int numberOfOutputs,
													const std::vector<unsigned int> & userGivenTimeStepList);

public:
	// Method to read the input parameters from a file and load them into the class member variables
	void loadInputParameters(dealii::ParameterHandler & parameter_handler, unsigned int _number_of_variables,
								unsigned int _number_of_materials, unsigned int _number_of_pp_variables,
								unsigned int _number_of_constants);

	// Method to create the list of BCs from the user input strings (called from loadInputParameters)
	void load_BC_list(std::vector<std::string> list_of_BCs, std::vector<varBCs<dim> > & BC_list);



	// Meshing parameters
	std::vector<double> domain_size;
	std::vector<unsigned int> subdivisions;
	unsigned int refine_factor;
	unsigned int degree;

	// Mesh refinement parameters
	bool h_adaptivity;
	unsigned int max_refinement_level;
	unsigned int min_refinement_level;

	std::vector<int> refine_criterion_fields;
	std::vector<double> refine_window_max;
	std::vector<double> refine_window_min;

	unsigned int skip_remeshing_steps;

	// Output parameters (delete all but the last)
	unsigned int skip_print_steps;
	std::string output_file_type;
	std::string output_file_name;
	std::vector<unsigned int> outputTimeStepList;

	bool calc_energy;

	// Time step parameters
	double dtValue;
	double finalTime;
	unsigned int totalIncrements;

	// Elliptic solver parameters
	std::string solver_type;
	bool abs_tol;
	double solver_tolerance;
	unsigned int max_solver_iterations;

	// Variable inputs
	unsigned int number_of_variables;

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


	// Variables needed to calculate the LHS
	std::vector<variable_info> varInfoListRHS;
	std::vector<variable_info> resInfoListRHS;

	// Variables needed to calculate the LHS
	unsigned int num_var_LHS;
	std::vector<variable_info> varInfoListLHS;

	// Variables for loading in initial conditions
	std::vector<bool> load_ICs;
	std::vector<bool> load_serial_file;
	std::vector<std::string> load_file_name;
	std::vector<std::string> load_field_name;

	// Elasticity tensor
	std::vector<dealii::Tensor<2, 2*dim-1+dim/3, dealii::VectorizedArray<double> > > CIJ_list;
	list_of_CIJ<dim> material_moduli;

	bool nucleation_occurs;

	// Postprocessing parameters
	unsigned int pp_number_of_variables;
	bool postProcessingRequired;

	std::vector<std::string> pp_var_name;
	std::vector<std::string> pp_var_type;
	std::vector<std::string> pp_var_eq_type;

	std::vector<bool> pp_need_value;
	std::vector<bool> pp_need_gradient;
	std::vector<bool> pp_need_hessian;
	std::vector<bool> pp_value_residual;
	std::vector<bool> pp_gradient_residual;

	// Variable and residual info
	std::vector<variable_info> pp_varInfoList;
	std::vector<variable_info> pp_resInfoList;

	// List of boundary conditions
	std::vector<varBCs<dim> > BC_list;

	// List of user-defined constants
	std::vector<boost::variant<double,int,bool> > model_constants;

};

#endif /* INCLUDE_USERINPUTPARAMETERS_H_ */
