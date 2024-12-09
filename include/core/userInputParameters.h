// Class to load in the user input from parameters.h and the variable definition
// part of equations.h

#ifndef INCLUDE_USERINPUTPARAMETERS_H_
#define INCLUDE_USERINPUTPARAMETERS_H_

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/vector.h>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/unordered_map.hpp>
#include <boost/variant.hpp>

#include <core/boundary_conditions/varBCs.h>
#include <core/inputFileReader.h>
#include <core/model_variables.h>
#include <core/refinement/RefinementCriterion.h>
#include <core/solvers/SolverParameters.h>
#include <core/variableAttributeLoader.h>
#include <nucleation/nucleationParameters.h>
#include <unordered_map>
#include <vector>

enum elasticityModel
{
  ISOTROPIC,
  TRANSVERSE,
  ORTHOTROPIC,
  ANISOTROPIC,
  ANISOTROPIC2D
};

template <int dim>
class userInputParameters
{
public:
  /**
   * \brief Constructor. Reads in user input parameters from file and loads them into
   * member variables.
   */
  userInputParameters(inputFileReader          &input_file_reader,
                      dealii::ParameterHandler &parameter_handler,
                      variableAttributeLoader   variable_attributes);

  /**
   * \brief Creates a list of BCs to store in BC_list object.
   */
  void
  load_BC_list(const std::vector<std::string> &list_of_BCs);

  /**
   * \brief Assign the boundary condition to the varBC<dim> object given some boundary
   * condition string vector.
   */
  void
  assign_boundary_conditions(std::vector<std::string> &boundary_condition_list,
                             varBCs<dim>              &boundary_condition);

  // Map linking the model constant name to its index
  std::unordered_map<std::string, unsigned int> model_constant_name_map;

  /**
   * \brief Retrieve the double from the `model_constants` that are defined from the
   * parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * \param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] double
  get_model_constant_double(const std::string &constant_name) const
  {
    Assert(model_constant_name_map.find(constant_name) != model_constant_name_map.end(),
           dealii::ExcMessage(
             "PRISMS-PF Error: Mismatch between constants in parameters.prm and "
             "customPDE.h. The constant that you attempted to access was " +
             constant_name + "."));

    return boost::get<double>(model_constants[model_constant_name_map.at(constant_name)]);
  };

  /**
   * \brief Retrieve the int from the `model_constants` that are defined from the
   * parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * \param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] int
  get_model_constant_int(const std::string &constant_name) const
  {
    Assert(model_constant_name_map.find(constant_name) != model_constant_name_map.end(),
           dealii::ExcMessage(
             "PRISMS-PF Error: Mismatch between constants in parameters.prm and "
             "customPDE.h. The constant that you attempted to access was " +
             constant_name + "."));

    return boost::get<int>(model_constants[model_constant_name_map.at(constant_name)]);
  };

  /**
   * \brief Retrieve the bool from the `model_constants` that are defined from the
   * parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * \param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] bool
  get_model_constant_bool(const std::string &constant_name) const
  {
    Assert(model_constant_name_map.find(constant_name) != model_constant_name_map.end(),
           dealii::ExcMessage(
             "PRISMS-PF Error: Mismatch between constants in parameters.prm and "
             "customPDE.h. The constant that you attempted to access was " +
             constant_name + "."));

    return boost::get<bool>(model_constants[model_constant_name_map.at(constant_name)]);
  };

  /**
   * \brief Retrieve the rank 1 tensor from the `model_constants` that are defined from
   * the parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * \param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] dealii::Tensor<1, dim>
  get_model_constant_rank_1_tensor(const std::string &constant_name) const
  {
    Assert(model_constant_name_map.find(constant_name) != model_constant_name_map.end(),
           dealii::ExcMessage(
             "PRISMS-PF Error: Mismatch between constants in parameters.prm and "
             "customPDE.h. The constant that you attempted to access was " +
             constant_name + "."));

    return boost::get<dealii::Tensor<1, dim>>(
      model_constants[model_constant_name_map.at(constant_name)]);
  };

  /**
   * \brief Retrieve the rank 2 tensor from the `model_constants` that are defined from
   * the parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * \param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] dealii::Tensor<2, dim>
  get_model_constant_rank_2_tensor(const std::string &constant_name) const
  {
    Assert(model_constant_name_map.find(constant_name) != model_constant_name_map.end(),
           dealii::ExcMessage(
             "PRISMS-PF Error: Mismatch between constants in parameters.prm and "
             "customPDE.h. The constant that you attempted to access was " +
             constant_name + "."));

    return boost::get<dealii::Tensor<2, dim>>(
      model_constants[model_constant_name_map.at(constant_name)]);
  };

  /**
   * \brief Retrieve the elasticity tensor from the `model_constants` that are defined
   * from the parameters.prm parser. This is essentially just a wrapper for boost::get.
   *
   * \param constant_name Name of the constant to retrieve.
   */
  [[nodiscard]] dealii::Tensor<2, 2 * dim - 1 + dim / 3>
  get_model_constant_elasticity_tensor(const std::string &constant_name) const
  {
    Assert(model_constant_name_map.find(constant_name) != model_constant_name_map.end(),
           dealii::ExcMessage(
             "PRISMS-PF Error: Mismatch between constants in parameters.prm and "
             "customPDE.h. The constant that you attempted to access was " +
             constant_name + "."));

    return boost::get<dealii::Tensor<2, 2 * dim - 1 + dim / 3>>(
      model_constants[model_constant_name_map.at(constant_name)]);
  };

  // Method to load in the variable attributes
  void
  loadVariableAttributes(const variableAttributeLoader &variable_attributes);

  // Nucleation attribute methods
  [[nodiscard]] std::vector<double>
  get_nucleus_semiaxes(unsigned int var_index) const
  {
    return nucleation_parameters_list[nucleation_parameters_list_index.at(var_index)]
      .semiaxes;
  };

  [[nodiscard]] std::vector<double>
  get_nucleus_freeze_semiaxes(unsigned int var_index) const
  {
    return nucleation_parameters_list[nucleation_parameters_list_index.at(var_index)]
      .freeze_semiaxes;
  };

  [[nodiscard]] std::vector<double>
  get_nucleus_rotation(unsigned int var_index) const
  {
    return nucleation_parameters_list[nucleation_parameters_list_index.at(var_index)]
      .ellipsoid_rotation;
  };

  [[nodiscard]] double
  get_no_nucleation_border_thickness(unsigned int var_index) const
  {
    return nucleation_parameters_list[nucleation_parameters_list_index.at(var_index)]
      .no_nucleation_border_thickness;
  };

  [[nodiscard]] double
  get_nucleus_hold_time(unsigned int var_index) const
  {
    return nucleation_parameters_list[nucleation_parameters_list_index.at(var_index)]
      .hold_time;
  };

  dealii::Tensor<2, dim, double>
  get_nucleus_rotation_matrix(unsigned int var_index) const
  {
    return nucleation_parameters_list[nucleation_parameters_list_index.at(var_index)]
      .rotation_matrix;
  };

  // Meshing parameters
  std::vector<double>       domain_size;
  std::vector<unsigned int> subdivisions;
  unsigned int              refine_factor;
  unsigned int              degree;

  // Mesh refinement parameters
  bool                             h_adaptivity;
  unsigned int                     max_refinement_level;
  unsigned int                     min_refinement_level;
  unsigned int                     skip_remeshing_steps;
  std::vector<RefinementCriterion> refinement_criteria;

  // Output parameters
  unsigned int              skip_print_steps;
  std::string               output_file_type;
  bool                      output_vtu_per_process;
  std::string               output_file_name;
  std::vector<unsigned int> outputTimeStepList;
  bool                      print_timing_with_output;

  // Time step parameters
  double       dtValue;
  double       finalTime;
  unsigned int totalIncrements;

  // Elliptic solver parameters
  LinearSolverParameters linear_solver_parameters;

  // Nonlinear solver parameters
  NonlinearSolverParameters nonlinear_solver_parameters;

  // Pinning point parameters
  boost::unordered_map<unsigned int, dealii::Point<dim>> pinned_point;

  // Variable inputs (I might be able to leave some/all of these in
  // variable_attributes)
  unsigned int number_of_variables;

  // Variables needed to calculate the RHS
  unsigned int               num_var_explicit_RHS, num_var_nonexplicit_RHS;
  std::vector<variable_info> varInfoListExplicitRHS, varInfoListNonexplicitRHS;

  // Variables needed to calculate the LHS
  unsigned int               num_var_LHS;
  std::vector<variable_info> varInfoListLHS;
  std::vector<variable_info> varChangeInfoListLHS;

  // Variables for loading in initial conditions
  std::vector<bool>        load_ICs;
  std::vector<bool>        load_parallel_file;
  std::vector<std::string> load_file_name;
  std::vector<std::string> load_field_name;

  // Variables for saving/loading checkpoints
  bool                      resume_from_checkpoint;
  std::vector<unsigned int> checkpointTimeStepList;

  // Postprocessing parameters
  unsigned int              pp_number_of_variables;
  unsigned int              num_integrated_fields;
  bool                      postProcessingRequired;
  std::vector<unsigned int> integrated_field_indices;

  // Variable and residual info
  std::vector<variable_info> pp_varInfoList;
  std::vector<variable_info> pp_baseVarInfoList;

  // List of boundary conditions
  std::vector<varBCs<dim>> BC_list;

  // List of user-defined constants
  std::vector<boost::variant<double,
                             int,
                             bool,
                             dealii::Tensor<1, dim>,
                             dealii::Tensor<2, dim>,
                             dealii::Tensor<2, 2 * dim - 1 + dim / 3>>>
    model_constants;

  // Nucleation parameters
  bool                      nucleation_occurs;
  std::vector<unsigned int> nucleating_variable_indices;
  std::vector<unsigned int> nucleation_need_value;
  bool                      evolution_before_nucleation;
  // Declare later
  // bool multiple_nuclei_per_order_parameter;
  double min_distance_between_nuclei; // Only enforced for nuclei placed during
                                      // the same time step
  double       nucleation_order_parameter_cutoff;
  unsigned int steps_between_nucleation_attempts;
  double       nucleation_start_time;
  double       nucleation_end_time;

  // Grain remapping parameters
  bool                      grain_remapping_activated;
  std::vector<unsigned int> variables_for_remapping; // Note: this should be a sorted list
  unsigned int              skip_grain_reassignment_steps;
  double                    order_parameter_threshold;
  double                    buffer_between_grains;

  bool        load_grain_structure;
  std::string load_vtk_file_type; // adding this string to know what type of vtk file you
                                  // want to read, it will be passed to
                                  // initialconditions.cc
  double       min_radius_for_loading_grains;
  std::string  grain_structure_filename;
  std::string  grain_structure_variable_name;
  unsigned int num_grain_smoothing_cycles;

private:
  /**
   * \brief Assign the provided user inputs to parameters for anything related to the
   * spatial discretiziation.
   */
  void
  assign_spatial_discretization_parameters(dealii::ParameterHandler &parameter_handler,
                                           variableAttributeLoader  &variable_attributes);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to the
   * temporal discretiziation.
   */
  void
  assign_temporal_discretization_parameters(dealii::ParameterHandler &parameter_handler,
                                            variableAttributeLoader &variable_attributes);
  /**
   * \brief Assign the provided user inputs to parameters for anything related to linear
   * solves.
   */
  void
  assign_linear_solve_parameters(dealii::ParameterHandler &parameter_handler,
                                 variableAttributeLoader  &variable_attributes);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * nonlinear solves.
   */
  void
  assign_nonlinear_solve_parameters(dealii::ParameterHandler &parameter_handler,
                                    variableAttributeLoader  &variable_attributes);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * outputs.
   */
  void
  assign_output_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * loading in initial condition.
   */
  void
  assign_load_initial_condition_parameters(dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * nucleation.
   */
  void
  assign_nucleation_parameters(dealii::ParameterHandler &parameter_handler,
                               variableAttributeLoader  &variable_attributes);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * grain remapping and grain vtk load-in.
   */
  void
  assign_grain_parameters(dealii::ParameterHandler &parameter_handler,
                          variableAttributeLoader  &variable_attributes);

  /**
   * \brief Assign the provided user inputs to parameters for anything related to
   * boundary conditions.
   */
  void
  assign_boundary_condition_parameters(dealii::ParameterHandler &parameter_handler,
                                       variableAttributeLoader  &variable_attributes);

  // Method to create the list of time steps where the results should be output
  // (called from loadInputParameters)
  std::vector<unsigned int>
  setTimeStepList(const std::string               &outputSpacingType,
                  unsigned int                     numberOfOutputs,
                  const std::vector<unsigned int> &userGivenTimeStepList);

  void
  load_user_constants(inputFileReader          &input_file_reader,
                      dealii::ParameterHandler &parameter_handler);

  /**
   * \brief Compute the number of tensor rows.
   */
  unsigned int
  compute_tensor_parentheses(const unsigned int              n_elements,
                             const std::vector<std::string> &tensor_elements);

  /**
   * \brief Remove and leading and trailing parentheses.
   */
  void
  remove_parentheses(std::vector<std::string> &tensor_elements);

  /**
   * \brief Compute a 1st rank tensor from user inputs .
   */
  dealii::Tensor<1, dim>
  compute_rank_1_tensor_constant(const unsigned int       n_elements,
                                 std::vector<std::string> tensor_elements);

  /**
   * \brief Compute a 2nd rank tensor from user inputs .
   */
  dealii::Tensor<2, dim>
  compute_rank_2_tensor_constant(const unsigned int       n_elements,
                                 std::vector<std::string> tensor_elements);

  /**
   * \brief Assign the specified user constant to whatever type.
   */
  void
  assign_user_constant(std::vector<std::string> &model_constants_strings);

  /**
   * \brief Assign the primitive user constants (e.g., int, double, bool).
   */
  void
  assign_primitive_user_constant(std::vector<std::string> &model_constants_strings);

  [[nodiscard]] dealii::Tensor<2, 2 * dim - 1 + dim / 3>
  get_Cij_tensor(std::vector<double> elastic_constants,
                 const std::string  &elastic_const_symmetry) const;

  dealii::Tensor<2, 2 * dim - 1 + dim / 3>
  getCIJMatrix(const elasticityModel       model,
               const std::vector<double>  &constants,
               dealii::ConditionalOStream &pcout) const;

  // Private nucleation variables
  std::vector<nucleationParameters<dim>> nucleation_parameters_list;
  std::map<unsigned int, unsigned int>   nucleation_parameters_list_index;
};

#endif /* INCLUDE_USERINPUTPARAMETERS_H_ */
