// Class to load in the user input from parameters.h and the variable definition
// part of equations.h

#ifndef INCLUDE_LOADINPUTS_H_
#define INCLUDE_LOADINPUTS_H_

class userInputParameters
{
public:
	// Method to read user input from parameters.h or equations.h and load it into the class member variables
	void loadUserInput();

	// Input parameters
	const static unsigned int dim = problemDIM;

	// Meshing parameters
	std::vector<double> domain_size;
	std::vector<unsigned int> subdivisions;
	unsigned int refine_factor;

	// Mesh refinement parameters
	bool h_adaptivity;
	unsigned int max_refinement_level;
	unsigned int min_refinement_level;

	std::vector<int> refine_criterion_fields;
	std::vector<double> refine_window_max;
	std::vector<double> refine_window_min;

	unsigned int skip_remeshing_steps;

	// Output parameters
	bool write_output;
	std::string output_condition;
	unsigned int num_outputs;
	std::vector<unsigned int> user_given_time_step_list;
	unsigned int skip_print_steps;
	std::string output_file_type;

	bool calc_energy;

	// Time step parameters
	double dtValue;
	double finalTime;
	unsigned int totalIncrements;

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

	bool nucleation_occurs;

};

#endif /* INCLUDE_LOADINPUTS_H_ */
