// Class to load in the user input from parameters.h and the variable definition
// part of equations.h

#ifndef INCLUDE_USERINPUTPARAMETERS_H_
#define INCLUDE_USERINPUTPARAMETERS_H_


#include <deal.II/lac/vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include "model_variables.h"
#include "varBCs.h"
#include "RefinementCriterion.h"
#include "inputFileReader.h"
#include "varTypeEnums.h"
#include "variableAttributeLoader.h"
#include "nucleationParameters.h"
#include "SolverParameters.h"
#include <deal.II/base/conditional_ostream.h>
#include <boost/variant.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <vector>
#include <iostream>
#include <unordered_map>

enum elasticityModel {ISOTROPIC, TRANSVERSE, ORTHOTROPIC, ANISOTROPIC, ANISOTROPIC2D};

template <int dim>
class userInputParameters
{

public:
	// Method to read the input parameters from a file and load them into the class member variables
	userInputParameters(inputFileReader & input_file_reader, dealii::ParameterHandler & parameter_handler, variableAttributeLoader variable_attributes);

	// Method to create the list of BCs from the user input strings (called from the constructor)
	void load_BC_list(const std::vector<std::string> list_of_BCs);

	// Map linking the model constant name to its index
	std::unordered_map<std::string,unsigned int> model_constant_name_map;

	// Methods to access members of 'model_constant', one for each type (since one can't template based on return values)
	// These are really just wrappers for Boost's 'get' function
	double get_model_constant_double(const std::string constant_name) const {return boost::get<double>(model_constants[model_constant_name_map.at(constant_name)]);};
	int get_model_constant_int(const std::string constant_name) const {return boost::get<int>(model_constants[model_constant_name_map.at(constant_name)]);};
	bool get_model_constant_bool(const std::string constant_name) const {return boost::get<bool>(model_constants[model_constant_name_map.at(constant_name)]);};
	dealii::Tensor<1,dim> get_model_constant_rank_1_tensor(const std::string constant_name) const {return boost::get<dealii::Tensor<1,dim> >(model_constants[model_constant_name_map.at(constant_name)]);};
	dealii::Tensor<2,dim> get_model_constant_rank_2_tensor(const std::string constant_name) const {return boost::get<dealii::Tensor<2,dim> >(model_constants[model_constant_name_map.at(constant_name)]);};
	dealii::Tensor<2,2*dim-1+dim/3> get_model_constant_elasticity_tensor(const std::string constant_name) const {return boost::get<dealii::Tensor<2,2*dim-1+dim/3> >(model_constants[model_constant_name_map.at(constant_name)]);};

	// Method to load in the variable attributes
	void loadVariableAttributes(variableAttributeLoader variable_attributes);

	// Nucleation attribute methods
	std::vector<double> get_nucleus_semiaxes(unsigned int var_index) const { return nucleation_parameters_list[nucleation_parameters_list_index.at(var_index)].semiaxes; };
	std::vector<double> get_nucleus_freeze_semiaxes(unsigned int var_index) const { return nucleation_parameters_list[nucleation_parameters_list_index.at(var_index)].freeze_semiaxes; };
	std::vector<double> get_nucleus_rotation(unsigned int var_index) const { return nucleation_parameters_list[nucleation_parameters_list_index.at(var_index)].ellipsoid_rotation; };
	double get_no_nucleation_border_thickness(unsigned int var_index) const { return nucleation_parameters_list[nucleation_parameters_list_index.at(var_index)].no_nucleation_border_thickness; };
	double get_nucleus_hold_time(unsigned int var_index) const { return nucleation_parameters_list[nucleation_parameters_list_index.at(var_index)].hold_time; };
	dealii::Tensor<2,dim,double> get_nucleus_rotation_matrix(unsigned int var_index) const { return nucleation_parameters_list[nucleation_parameters_list_index.at(var_index)].rotation_matrix; };

	// Meshing parameters
	std::vector<double> domain_size;
	std::vector<unsigned int> subdivisions;
	unsigned int refine_factor;
	unsigned int degree;

	// Mesh refinement parameters
	bool h_adaptivity;
	unsigned int max_refinement_level;
	unsigned int min_refinement_level;
    unsigned int skip_remeshing_steps;
    std::vector<RefinementCriterion> refinement_criteria;

	// Output parameters
	unsigned int skip_print_steps;
	std::string output_file_type;
	bool output_vtu_per_process;
	std::string output_file_name;
	std::vector<unsigned int> outputTimeStepList;
    bool print_timing_with_output;

	// Time step parameters
	double dtValue;
	double finalTime;
	unsigned int totalIncrements;

	// Elliptic solver parameters
    LinearSolverParameters linear_solver_parameters;

    // Nonlinear solver parameters
    NonlinearSolverParameters nonlinear_solver_parameters;

	// Variable inputs (I might be able to leave some/all of these in variable_attributes)
	unsigned int number_of_variables;

	std::vector<std::string> var_name;
	std::vector<fieldType> var_type;
	std::vector<PDEType> var_eq_type;

    std::vector<bool> var_nonlinear;

	// Variables needed to calculate the RHS
    unsigned int num_var_explicit_RHS, num_var_nonexplicit_RHS;
	std::vector<variable_info> varInfoListExplicitRHS, varInfoListNonexplicitRHS;

	// Variables needed to calculate the LHS
	unsigned int num_var_LHS;
	std::vector<variable_info> varInfoListLHS;
    std::vector<variable_info> varChangeInfoListLHS;

	// Variables for loading in initial conditions
	std::vector<bool> load_ICs;
	std::vector<bool> load_parallel_file;
	std::vector<std::string> load_file_name;
	std::vector<std::string> load_field_name;

	// Variables for saving/loading checkpoints
	bool resume_from_checkpoint;
	std::vector<unsigned int> checkpointTimeStepList;

	// Postprocessing parameters
	unsigned int pp_number_of_variables;
	unsigned int num_integrated_fields;
	bool postProcessingRequired;
	std::vector<bool> pp_calc_integral;
	std::vector<unsigned int> integrated_field_indices;

	std::vector<std::string> pp_var_name;
	std::vector<fieldType> pp_var_type;

	// Variable and residual info
	std::vector<variable_info> pp_varInfoList;
	std::vector<variable_info> pp_baseVarInfoList;

	// List of boundary conditions
	std::vector<varBCs<dim> > BC_list;

	// List of user-defined constants
	std::vector<boost::variant<double, int, bool,dealii::Tensor<1,dim>, dealii::Tensor<2,dim>, dealii::Tensor<2,2*dim-1+dim/3> > > model_constants;

	// Nucleation parameters
	bool nucleation_occurs;
	std::vector<unsigned int> nucleating_variable_indices;
	std::vector<unsigned int> nucleation_need_value;

  bool multiple_nuclei_per_order_parameter;
	double min_distance_between_nuclei; // Only enforced for nuclei placed during the same time step
	double min_distance_between_OP;
	double nucleation_order_parameter_cutoff;
	unsigned int steps_between_nucleation_attempts;
  double nucleation_start_time;
  double nucleation_end_time;
  
    // Grain remapping parameters
    bool grain_remapping_activated;
    std::vector<unsigned int> variables_for_remapping; // Note: this should be a sorted list
    unsigned int skip_grain_reassignment_steps;
    double order_parameter_threshold;
    double buffer_between_grains;

    bool load_grain_structure;
    bool load_unstructured_grid;
    double min_radius_for_loading_grains;
    std::string grain_structure_filename;
    std::string grain_structure_variable_name;
    unsigned int num_grain_smoothing_cycles;

private:
	// Method to create the list of time steps where the results should be output (called from loadInputParameters)
	std::vector<unsigned int> setTimeStepList(const std::string outputSpacingType, unsigned int numberOfOutputs,
													const std::vector<unsigned int> & userGivenTimeStepList);

	void load_user_constants(inputFileReader & input_file_reader, dealii::ParameterHandler & parameter_handler);

	dealii::Tensor<2,2*dim-1+dim/3> get_Cij_tensor(std::vector<double> elastic_constants, const std::string elastic_const_symmetry) const;

	dealii::Tensor<2,2*dim-1+dim/3> getCIJMatrix(const elasticityModel model, const std::vector<double> constants, dealii::ConditionalOStream & pcout) const;

	// Private nucleation variables
	std::vector<nucleationParameters<dim> > nucleation_parameters_list;
	std::map<unsigned int, unsigned int> nucleation_parameters_list_index;

};

#endif /* INCLUDE_USERINPUTPARAMETERS_H_ */
