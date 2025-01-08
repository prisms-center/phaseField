#include <config.h>
#include <core/type_enums.h>
#include <core/user_inputs/user_input_parameters.h>
#include <iomanip>

template <int dim>
userInputParameters<dim>::userInputParameters(inputFileReader          &input_file_reader,
                                              dealii::ParameterHandler &parameter_handler)
  : var_attributes(input_file_reader.var_attributes)
  , pp_attributes(input_file_reader.pp_attributes)
{
  // Spatial discretization
  assign_spatial_discretization_parameters(parameter_handler);

  // Time stepping parameters
  assign_temporal_discretization_parameters(parameter_handler);

  // Output parameters
  assign_output_parameters(parameter_handler);

  // Boundary parameters
  assign_boundary_parameters(parameter_handler);

  // Linear solve parameters
  assign_linear_solve_parameters(parameter_handler);

  // Load the user-defined constants
  load_model_constants(input_file_reader, parameter_handler);

  // Print all the parameters to summary.log
  spatial_discretization.print_parameter_summary();
  temporal_discretization.print_parameter_summary();
  output_parameters.print_parameter_summary();
  boundary_parameters.print_parameter_summary();
  linear_solve_parameters.print_parameter_summary();
}

template <int dim>
void
userInputParameters<dim>::assign_spatial_discretization_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  std::vector<std::string> axis_labels = {"X", "Y", "Z"};
  for (uint i = 0; i < dim; ++i)
    {
      spatial_discretization.domain_size[i] =
        parameter_handler.get_double("Domain size " + axis_labels[i]);
      spatial_discretization.subdivisions[i] =
        parameter_handler.get_integer("Subdivisions " + axis_labels[i]);
    }

  spatial_discretization.refine_factor = parameter_handler.get_integer("Refine factor");

  spatial_discretization.degree = parameter_handler.get_integer("Element degree");

  // Adaptive meshing parameters
  spatial_discretization.has_adaptivity = parameter_handler.get_bool("Mesh adaptivity");

  AssertThrow(
    (!spatial_discretization.has_adaptivity || dim != 1),
    dealii::ExcMessage(
      "Adaptive meshing for the matrix-free method is not currently supported."));

  spatial_discretization.remeshing_frequency =
    parameter_handler.get_integer("Steps between remeshing operations");

  spatial_discretization.max_refinement_level =
    parameter_handler.get_integer("Max refinement level");
  spatial_discretization.min_refinement_level =
    parameter_handler.get_integer("Min refinement level");

  // Enforce that the initial refinement level must be between the max and min
  // level
  if (spatial_discretization.has_adaptivity &&
      ((spatial_discretization.refine_factor <
        spatial_discretization.min_refinement_level) ||
       (spatial_discretization.refine_factor >
        spatial_discretization.max_refinement_level)))
    {
      std::cerr << "PRISMS-PF Error: The initial refinement factor must be "
                   "between the maximum and minimum refinement levels when "
                   "adaptive meshing is enabled.\n";
      std::cerr << "Initial refinement level: " << spatial_discretization.refine_factor
                << " Maximum and minimum refinement levels: "
                << spatial_discretization.max_refinement_level << ", "
                << spatial_discretization.min_refinement_level << "\n";
      abort();
    }

  // The adaptivity criterion for each variable has its own subsection
  for (const auto &[index, variable] : var_attributes)
    {
      std::string subsection_text = "Refinement criterion: ";
      subsection_text.append(variable.name);

      parameter_handler.enter_subsection(subsection_text);
      {
        const std::string crit_type_string = parameter_handler.get("Criterion type");
        if (!crit_type_string.empty())
          {
            RefinementCriterion new_criterion;
            new_criterion.variable_index = index;
            new_criterion.variable_name  = variable.name;
            if (boost::iequals(crit_type_string, "VALUE"))
              {
                new_criterion.criterion_type = criterion_value;
                new_criterion.value_lower_bound =
                  parameter_handler.get_double("Value lower bound");
                new_criterion.value_upper_bound =
                  parameter_handler.get_double("Value upper bound");

                // Check to make sure that the upper bound is greater than or
                // equal to the lower bound
                if (new_criterion.value_upper_bound < new_criterion.value_lower_bound)
                  {
                    std::cerr << "PRISMS-PF Error: The upper bound for "
                                 "refinement for variable "
                              << new_criterion.variable_name
                              << " is less than the lower bound. Please "
                                 "correct this in the parameters file.\n";
                  }
              }
            else if (boost::iequals(crit_type_string, "GRADIENT"))
              {
                new_criterion.criterion_type = criterion_gradient;
                new_criterion.gradient_lower_bound =
                  parameter_handler.get_double("Gradient magnitude lower bound");
              }
            else if (boost::iequals(crit_type_string, "VALUE_AND_GRADIENT"))
              {
                new_criterion.criterion_type = criterion_value | criterion_gradient;
                new_criterion.value_lower_bound =
                  parameter_handler.get_double("Value lower bound");
                new_criterion.value_upper_bound =
                  parameter_handler.get_double("Value upper bound");
                new_criterion.gradient_lower_bound =
                  parameter_handler.get_double("Gradient magnitude lower bound");

                // Check to make sure that the upper bound is greater than or
                // equal to the lower bound
                if (new_criterion.value_upper_bound < new_criterion.value_lower_bound)
                  {
                    std::cerr << "PRISMS-PF Error: The upper bound for "
                                 "refinement for variable "
                              << new_criterion.variable_name
                              << " is less than the lower bound. Please "
                                 "correct this in the parameters file.\n";
                  }
              }
            else
              {
                std::cerr << "PRISMS-PF Error: The refinement criteria type "
                             "found in the parameters file, "
                          << crit_type_string
                          << ", is not an allowed type. The allowed types are "
                             "VALUE, GRADIENT, VALUE_AND_GRADIENT\n";
                abort();
              }
            spatial_discretization.refinement_criteria.push_back(new_criterion);
          }
      }
      parameter_handler.leave_subsection();
    }
}

template <int dim>
void
userInputParameters<dim>::assign_temporal_discretization_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  temporal_discretization.dt = parameter_handler.get_double("Time step");
  temporal_discretization.final_time =
    parameter_handler.get_double("Simulation end time");
  temporal_discretization.total_increments =
    static_cast<uint>(parameter_handler.get_integer("Number of time steps"));

  // If all of the variables are `TIME_INDEPENDENT` or `AUXILIARY`, then total_increments
  // should be 1 and final_time should be 0
  bool only_time_independent_pdes = true;
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.pde_type == PDEType::EXPLICIT_TIME_DEPENDENT ||
          variable.pde_type == PDEType::IMPLICIT_TIME_DEPENDENT)
        {
          only_time_independent_pdes = false;
          break;
        }
    }

  if (only_time_independent_pdes)
    {
      temporal_discretization.total_increments = 1;
      return;
    }

  AssertThrow(temporal_discretization.dt > 0.0,
              dealii::ExcMessage(
                "The timestep (dt) must be greater than zero for transient problems."));

  // Pick the maximum specified time since the default values are zero
  temporal_discretization.final_time =
    std::max(temporal_discretization.final_time,
             temporal_discretization.dt * temporal_discretization.total_increments);
  temporal_discretization.total_increments = static_cast<uint>(
    std::ceil(temporal_discretization.final_time / temporal_discretization.dt));
}

template <int dim>
void
userInputParameters<dim>::assign_output_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  output_parameters.output_condition = parameter_handler.get("Output condition");
  output_parameters.n_outputs =
    static_cast<uint>(parameter_handler.get_integer("Number of outputs"));
  output_parameters.user_output_list =
    dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(
      parameter_handler.get("List of time steps to output")));
  output_parameters.output_frequency = parameter_handler.get_integer("Skip print steps");
  output_parameters.output_file_type = parameter_handler.get("Output file type");
  output_parameters.output_file_name = parameter_handler.get("Output file name (base)");
  output_parameters.print_timing_with_output =
    parameter_handler.get_bool("Print timing information with output");

  // Determine the list of increments where output is neccessary
  output_parameters.compute_output_list(temporal_discretization);
}

template <int dim>
void
userInputParameters<dim>::assign_boundary_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  // Load the boundary condition variables into list of BCs (where each element
  // of the vector is one component of one variable)
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.field_type == SCALAR)
        {
          std::string bc_text = "Boundary condition for variable ";
          bc_text.append(variable.name);
          boundary_parameters.BC_list[index].emplace(0, parameter_handler.get(bc_text));
        }
      else
        {
          std::vector<std::string> axis_labels = {"x", "y", "z"};
          for (uint i = 0; i < dim; i++)
            {
              std::string bc_text = "Boundary condition for variable ";
              bc_text.append(variable.name + ", " + axis_labels[i] + " component");
              boundary_parameters.BC_list[index].emplace(i,
                                                         parameter_handler.get(bc_text));
            }
        }
    }

  // Compute the boundary conditions from the unfiltered string and throw errors, if
  // applicable.
  boundary_parameters.compute_boundary_conditions(var_attributes);

  for (const auto &[index, variable] : var_attributes)
    {
      std::string pinning_text = "Pinning point: ";
      pinning_text.append(variable.name);
      parameter_handler.enter_subsection(pinning_text);

      // Skip if the value is the default INT_MAX
      if (parameter_handler.get_double("value") == 2147483647)
        {
          parameter_handler.leave_subsection();
          continue;
        }

      // Otherwise, fill out point and value
      std::vector<std::string> axis_labels = {"x", "y", "z"};
      dealii::Point<dim>       point;
      for (uint i = 0; i < dim; ++i)
        {
          point[i] = parameter_handler.get_double(axis_labels[i]);
        }
      boundary_parameters.pinned_point_list
        .emplace(index, std::make_pair(parameter_handler.get_double("value"), point));

      parameter_handler.leave_subsection();
    }
}

template <int dim>
void
userInputParameters<dim>::assign_linear_solve_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.pde_type == TIME_INDEPENDENT ||
          variable.pde_type == IMPLICIT_TIME_DEPENDENT)
        {
          std::string subsection_text = "Linear solver parameters: ";
          subsection_text.append(variable.name);
          parameter_handler.enter_subsection(subsection_text);

          // Set the tolerance type
          const std::string type_string = parameter_handler.get("Tolerance type");
          if (boost::iequals(type_string, "ABSOLUTE_RESIDUAL"))
            {
              linear_solve_parameters.linear_solve[index].tolerance_type =
                solverToleranceType::ABSOLUTE_RESIDUAL;
            }
          else if (boost::iequals(type_string, "RELATIVE_RESIDUAL_CHANGE"))
            {
              linear_solve_parameters.linear_solve[index].tolerance_type =
                solverToleranceType::RELATIVE_RESIDUAL_CHANGE;
            }

          // Set the tolerance value
          linear_solve_parameters.linear_solve[index].tolerance =
            parameter_handler.get_double("Tolerance value");

          // Set the maximum number of iterations
          linear_solve_parameters.linear_solve[index].max_iterations =
            parameter_handler.get_integer("Maximum linear solver iterations");

          // Set preconditioner type and related parameters
          linear_solve_parameters.linear_solve[index].preconditioner =
            boost::iequals(parameter_handler.get("Preconditioner"), "GMG")
              ? preconditionerType::GMG
              : preconditionerType::NONE;

          linear_solve_parameters.linear_solve[index].smoothing_range =
            parameter_handler.get_double("Smoothing range");

          linear_solve_parameters.linear_solve[index].smoother_iterations =
            parameter_handler.get_integer("Smoother iterations");

          linear_solve_parameters.linear_solve[index].eig_cg_n_iterations =
            parameter_handler.get_integer("Eigenvalue CG iterations");

          parameter_handler.leave_subsection();
        }
    }
}

template <int dim>
void
userInputParameters<dim>::load_model_constants(
  inputFileReader          &input_file_reader,
  dealii::ParameterHandler &parameter_handler)
{
  for (const std::string &constant_name : input_file_reader.model_constant_names)
    {
      std::string constants_text = "Model constant ";
      constants_text.append(constant_name);

      std::vector<std::string> model_constants_strings =
        dealii::Utilities::split_string_list(parameter_handler.get(constants_text));

      user_constants.model_constants[constant_name] =
        user_constants.construct_user_constant(model_constants_strings);
    }
}

INSTANTIATE_UNI_TEMPLATE(userInputParameters)