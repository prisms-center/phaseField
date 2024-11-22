// Methods for the userInputParameters class
#include "../../include/userInputParameters.h"

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <cstddef>
#include <vector>

template <int dim>
void
userInputParameters<dim>::assign_spatial_discretization_parameters(
  dealii::ParameterHandler &parameter_handler,
  variableAttributeLoader  &variable_attributes)
{
  // Domain size & subdivisions
  domain_size.push_back(parameter_handler.get_double("Domain size X"));
  subdivisions.push_back(parameter_handler.get_integer("Subdivisions X"));

  if (dim > 1)
    {
      domain_size.push_back(parameter_handler.get_double("Domain size Y"));
      subdivisions.push_back(parameter_handler.get_integer("Subdivisions Y"));

      if (dim > 2)
        {
          domain_size.push_back(parameter_handler.get_double("Domain size Z"));
          subdivisions.push_back(parameter_handler.get_integer("Subdivisions Z"));
        }
    }

  refine_factor = parameter_handler.get_integer("Refine factor");

  degree = parameter_handler.get_integer("Element degree");

  // Adaptive meshing parameters
  h_adaptivity = parameter_handler.get_bool("Mesh adaptivity");
  skip_remeshing_steps =
    parameter_handler.get_integer("Steps between remeshing operations");

  max_refinement_level = parameter_handler.get_integer("Max refinement level");
  min_refinement_level = parameter_handler.get_integer("Min refinement level");

  // Enforce that the initial refinement level must be between the max and min
  // level
  if (h_adaptivity &&
      ((refine_factor < min_refinement_level) || (refine_factor > max_refinement_level)))
    {
      std::cerr << "PRISMS-PF Error: The initial refinement factor must be "
                   "between the maximum and minimum refinement levels when "
                   "adaptive meshing is enabled.\n";
      std::cerr << "Initial refinement level: " << refine_factor
                << " Maximum and minimum refinement levels: " << max_refinement_level
                << ", " << min_refinement_level << "\n";
      abort();
    }

  // The adaptivity criterion for each variable has its own subsection
  for (const auto &[index, variable] : variable_attributes.attributes)
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
            refinement_criteria.push_back(new_criterion);
          }
      }
      parameter_handler.leave_subsection();
    }
}

template <int dim>
void
userInputParameters<dim>::assign_temporal_discretization_parameters(
  dealii::ParameterHandler &parameter_handler,
  variableAttributeLoader  &variable_attributes)
{
  dtValue = parameter_handler.get_double("Time step");
  const int totalIncrements_temp =
    static_cast<int>(parameter_handler.get_integer("Number of time steps"));
  finalTime = parameter_handler.get_double("Simulation end time");

  // If all of the variables are ELLIPTIC, then totalIncrements should be 1 and
  // finalTime should be 0
  bool only_time_independent_pdes = true;
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (variable.eq_type == EXPLICIT_TIME_DEPENDENT ||
          variable.eq_type == IMPLICIT_TIME_DEPENDENT)
        {
          only_time_independent_pdes = false;
          break;
        }
    }

  // Determine the maximum number of time steps
  if (only_time_independent_pdes)
    {
      totalIncrements = 1;
      finalTime       = 0.0;
    }
  else
    {
      if ((totalIncrements_temp >= 0) && (finalTime >= 0.0))
        {
          if (std::ceil(finalTime / dtValue) > totalIncrements_temp)
            {
              totalIncrements = totalIncrements_temp;
              finalTime       = totalIncrements * dtValue;
            }
          else
            {
              totalIncrements = std::ceil(finalTime / dtValue);
            }
        }
      else if ((totalIncrements_temp >= 0) && (finalTime < 0.0))
        {
          totalIncrements = totalIncrements_temp;
          finalTime       = totalIncrements * dtValue;
        }
      else if ((totalIncrements_temp < 0) && (finalTime >= 0.0))
        {
          totalIncrements = std::ceil(finalTime / dtValue);
        }
      else
        {
          // Should change to an exception
          std::cerr << "Invalid selections for the final time and the number "
                       "of increments. At least one should be given in the "
                       "input file and should be positive.\n";
          std::cout << finalTime << " " << totalIncrements_temp << "\n";
          abort();
        }
    }
}

template <int dim>
void
userInputParameters<dim>::assign_linear_solve_parameters(
  dealii::ParameterHandler &parameter_handler,
  variableAttributeLoader  &variable_attributes)
{
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (variable.eq_type == TIME_INDEPENDENT ||
          variable.eq_type == IMPLICIT_TIME_DEPENDENT)
        {
          std::string subsection_text = "Linear solver parameters: ";
          subsection_text.append(variable.name);

          parameter_handler.enter_subsection(subsection_text);
          {
            // Set the tolerance type
            SolverToleranceType temp_type   = ABSOLUTE_RESIDUAL;
            const std::string   type_string = parameter_handler.get("Tolerance type");
            if (boost::iequals(type_string, "ABSOLUTE_RESIDUAL"))
              {
                temp_type = ABSOLUTE_RESIDUAL;
              }
            else if (boost::iequals(type_string, "RELATIVE_RESIDUAL_CHANGE"))
              {
                temp_type = RELATIVE_RESIDUAL_CHANGE;
              }
            else if (boost::iequals(type_string, "ABSOLUTE_SOLUTION_CHANGE"))
              {
                temp_type = ABSOLUTE_SOLUTION_CHANGE;
                std::cerr << "PRISMS-PF Error: Linear solver tolerance type "
                          << type_string
                          << " is not currently implemented, please use either "
                             "ABSOLUTE_RESIDUAL or RELATIVE_RESIDUAL_CHANGE\n";
                abort();
              }
            else
              {
                std::cerr << "PRISMS-PF Error: Linear solver tolerance type "
                          << type_string
                          << " is not one of the allowed values (ABSOLUTE_RESIDUAL, "
                             "RELATIVE_RESIDUAL_CHANGE, ABSOLUTE_SOLUTION_CHANGE)\n";
                abort();
              }

            // Set the tolerance value
            const double temp_value = parameter_handler.get_double("Tolerance value");

            // Set the maximum number of iterations
            const unsigned int temp_max_iterations =
              parameter_handler.get_integer("Maximum linear solver iterations");

            linear_solver_parameters.loadParameters(index,
                                                    temp_type,
                                                    temp_value,
                                                    temp_max_iterations);
          }
          parameter_handler.leave_subsection();
        }
    }
}

template <int dim>
void
userInputParameters<dim>::assign_nonlinear_solve_parameters(
  dealii::ParameterHandler &parameter_handler,
  variableAttributeLoader  &variable_attributes)
{
  nonlinear_solver_parameters.setMaxIterations(
    parameter_handler.get_integer("Maximum nonlinear solver iterations"));

  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (variable.is_nonlinear)
        {
          std::string subsection_text = "Nonlinear solver parameters: ";
          subsection_text.append(variable.name);

          parameter_handler.enter_subsection(subsection_text);
          {
            // Set the tolerance type
            SolverToleranceType temp_type   = ABSOLUTE_RESIDUAL;
            const std::string   type_string = parameter_handler.get("Tolerance type");
            if (boost::iequals(type_string, "ABSOLUTE_RESIDUAL"))
              {
                temp_type = ABSOLUTE_RESIDUAL;
              }
            else if (boost::iequals(type_string, "RELATIVE_RESIDUAL_CHANGE"))
              {
                temp_type = RELATIVE_RESIDUAL_CHANGE;
              }
            else if (boost::iequals(type_string, "ABSOLUTE_SOLUTION_CHANGE"))
              {
                temp_type = ABSOLUTE_SOLUTION_CHANGE;
              }
            else
              {
                std::cerr << "PRISMS-PF Error: Nonlinear solver tolerance type "
                          << type_string
                          << " is not one of the allowed values (ABSOLUTE_RESIDUAL, "
                             "RELATIVE_RESIDUAL_CHANGE, ABSOLUTE_SOLUTION_CHANGE)\n";
                abort();
              }

            // Set the tolerance value
            const double temp_value = parameter_handler.get_double("Tolerance value");

            // Set the backtrace damping flag
            const bool temp_backtrack_damping =
              parameter_handler.get_bool("Use backtracking line search damping");

            // Set the backtracking step size modifier
            const double temp_step_modifier =
              parameter_handler.get_double("Backtracking step size modifier");

            // Set the constant that determines how much the residual must
            // decrease to be accepted as sufficient
            const double temp_residual_decrease_coeff =
              parameter_handler.get_double("Backtracking residual decrease coefficient");

            // Set the default damping coefficient (used if backtracking isn't
            // used)
            const double temp_damping_coefficient =
              parameter_handler.get_double("Constant damping value");

            // Set whether to use the solution of Laplace's equation instead of
            // the IC in ICs_and_BCs.h as the initial guess for nonlinear, time
            // independent equations
            bool temp_laplace_for_initial_guess = false;
            if (variable.eq_type == TIME_INDEPENDENT)
              {
                temp_laplace_for_initial_guess = parameter_handler.get_bool(
                  "Use Laplace's equation to determine the initial guess");
              }
            else
              {
                if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
                  {
                    std::cout << "PRISMS-PF Warning: Laplace's equation is only used "
                                 "to generate the initial guess for time independent "
                                 "equations. The equation for variable "
                              << variable.name
                              << " is not a time independent equation. No initial "
                                 "guess is needed for this equation.\n";
                  }
              }

            nonlinear_solver_parameters.loadParameters(index,
                                                       temp_type,
                                                       temp_value,
                                                       temp_backtrack_damping,
                                                       temp_step_modifier,
                                                       temp_residual_decrease_coeff,
                                                       temp_damping_coefficient,
                                                       temp_laplace_for_initial_guess);
          }
          parameter_handler.leave_subsection();
        }
    }

  // Set the max number of nonlinear iterations
  bool any_nonlinear = false;
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      any_nonlinear |= variable.is_nonlinear;
    }
  if (!any_nonlinear)
    {
      nonlinear_solver_parameters.setMaxIterations(0);
    }
}

template <int dim>
void
userInputParameters<dim>::assign_output_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  const std::string      output_condition = parameter_handler.get("Output condition");
  const unsigned int     num_outputs = parameter_handler.get_integer("Number of outputs");
  const std::vector<int> user_given_time_step_list_temp =
    dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(
      parameter_handler.get("List of time steps to output")));
  std::vector<unsigned int> user_given_time_step_list;
  user_given_time_step_list.reserve(user_given_time_step_list_temp.size());
  for (const auto &time_step : user_given_time_step_list_temp)
    {
      user_given_time_step_list.push_back(time_step);
    }

  skip_print_steps = parameter_handler.get_integer("Skip print steps");
  output_file_type = parameter_handler.get("Output file type");
  output_file_name = parameter_handler.get("Output file name (base)");

  output_vtu_per_process =
    parameter_handler.get_bool("Output separate files per process");
  if ((output_file_type == "vtk") && (!output_vtu_per_process))
    {
      output_vtu_per_process = true;
      if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        {
          std::cout << "PRISMS-PF Warning: 'Output file type' given as 'vtk' and "
                       "'Output separate files per process' given as 'false'. Shared "
                       "output files are not supported for the vtk output format. "
                       "Separate files per process will be created.\n";
        }
    }

  print_timing_with_output =
    parameter_handler.get_bool("Print timing information with output");

  // Field variable definitions

  // Use these inputs to create a list of time steps where the code should
  // output, stored in the member
  outputTimeStepList =
    setTimeStepList(output_condition, num_outputs, user_given_time_step_list);

  // Parameters for checkpoint/restart
  resume_from_checkpoint = parameter_handler.get_bool("Load from a checkpoint");
  const std::string  checkpoint_condition = parameter_handler.get("Checkpoint condition");
  const unsigned int num_checkpoints =
    parameter_handler.get_integer("Number of checkpoints");

  const std::vector<int> user_given_checkpoint_time_step_list_temp =
    dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(
      parameter_handler.get("List of time steps to save checkpoints")));
  std::vector<unsigned int> user_given_checkpoint_time_step_list;
  user_given_checkpoint_time_step_list.reserve(
    user_given_checkpoint_time_step_list_temp.size());
  for (const auto &checkpoint_step : user_given_checkpoint_time_step_list_temp)
    {
      user_given_checkpoint_time_step_list.push_back(checkpoint_step);
    }
  checkpointTimeStepList = setTimeStepList(checkpoint_condition,
                                           num_checkpoints,
                                           user_given_checkpoint_time_step_list);
}

template <int dim>
void
userInputParameters<dim>::assign_load_initial_condition_parameters(
  dealii::ParameterHandler &parameter_handler)
{ // Variables for loading in PField ICs
  std::vector<std::string> load_ICs_temp = dealii::Utilities::split_string_list(
    parameter_handler.get("Load initial conditions"));
  std::vector<std::string> load_parallel_file_temp =
    dealii::Utilities::split_string_list(parameter_handler.get("Load parallel file"));

  if (boost::iequals(load_ICs_temp.at(0), "void"))
    {
      for (unsigned int var = 0; var < number_of_variables; var++)
        {
          load_ICs.push_back(false);
          load_parallel_file.push_back(false);
        }
    }
  else
    {
      for (unsigned int var = 0; var < number_of_variables; var++)
        {
          if (boost::iequals(load_ICs_temp.at(var), "true"))
            {
              load_ICs.push_back(true);
            }
          else
            {
              load_ICs.push_back(false);
            }
          if (boost::iequals(load_parallel_file_temp.at(var), "true"))
            {
              load_parallel_file.push_back(true);
            }
          else
            {
              load_parallel_file.push_back(false);
            }
        }
    }

  load_file_name =
    dealii::Utilities::split_string_list(parameter_handler.get("File names"));
  load_field_name = dealii::Utilities::split_string_list(
    parameter_handler.get("Variable names in the files"));
}

template <int dim>
void
userInputParameters<dim>::assign_nucleation_parameters(
  dealii::ParameterHandler &parameter_handler,
  variableAttributeLoader  &variable_attributes)
{
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (variable.nucleating_variable)
        {
          std::string nucleation_text = "Nucleation parameters: ";
          nucleation_text.append(variable.name);

          parameter_handler.enter_subsection(nucleation_text);
          {
            const unsigned int        var_index = index;
            const std::vector<double> semiaxes =
              dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(
                parameter_handler.get("Nucleus semiaxes (x, y, z)")));
            const std::vector<double> ellipsoid_rotation =
              dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(
                parameter_handler.get("Nucleus rotation in degrees (x, y, z)")));
            const std::vector<double> freeze_semiaxes =
              dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(
                parameter_handler.get("Freeze zone semiaxes (x, y, z)")));
            const double hold_time =
              parameter_handler.get_double("Freeze time following nucleation");
            const double no_nucleation_border_thickness =
              parameter_handler.get_double("Nucleation-free border thickness");

            const nucleationParameters<dim> temp(var_index,
                                                 semiaxes,
                                                 freeze_semiaxes,
                                                 ellipsoid_rotation,
                                                 hold_time,
                                                 no_nucleation_border_thickness);
            nucleation_parameters_list.push_back(temp);

            // Validate nucleation input
            if (semiaxes.size() < dim || semiaxes.size() > 3)
              {
                std::cerr << "PRISMS-PF Error: The number of nucleus semiaxes given in "
                             "the 'parameters.in' file must be at least the number of "
                             "dimensions and no more than 3.\n";
                abort();
              }
            if (freeze_semiaxes.size() < dim || freeze_semiaxes.size() > 3)
              {
                std::cerr << "PRISMS-PF Error: The number of nucleation freeze zone "
                             "semiaxes given in the 'parameters.in' file must be at "
                             "least the number of dimensions and no more than 3.\n";
                abort();
              }
            if (ellipsoid_rotation.size() != 3)
              {
                std::cerr << "PRISMS-PF Error: Exactly three nucleus rotation "
                             "angles must be given in the 'parameters.in' file.\n";
                abort();
              }
          }
          parameter_handler.leave_subsection();
        }
    }
  for (unsigned int i = 0; i < nucleation_parameters_list.size(); i++)
    {
      nucleation_parameters_list_index[nucleation_parameters_list.at(i).var_index] = i;
    }

  if (parameter_handler.get("Minimum allowed distance between nuclei") != "-1")
    {
      min_distance_between_nuclei =
        parameter_handler.get_double("Minimum allowed distance between nuclei");
    }
  else if (nucleation_parameters_list.size() > 1)
    {
      min_distance_between_nuclei =
        2.0 * (*(max_element(nucleation_parameters_list[0].semiaxes.begin(),
                             nucleation_parameters_list[0].semiaxes.end())));
    }
  evolution_before_nucleation =
    parameter_handler.get_bool("Enable evolution before nucleation");
  // Implement multiple order parameter nucleation later
  // multiple_nuclei_per_order_parameter = parameter_handler.get_bool("Allow
  // multiple nuclei per order parameter");
  nucleation_order_parameter_cutoff =
    parameter_handler.get_double("Order parameter cutoff value");
  steps_between_nucleation_attempts =
    parameter_handler.get_integer("Time steps between nucleation attempts");
  nucleation_start_time = parameter_handler.get_double("Nucleation start time");
  nucleation_end_time   = parameter_handler.get_double("Nucleation end time");
}

template <int dim>
void
userInputParameters<dim>::assign_grain_parameters(
  dealii::ParameterHandler &parameter_handler,
  variableAttributeLoader  &variable_attributes)
{ // Load the grain remapping parameters
  grain_remapping_activated = parameter_handler.get_bool("Activate grain reassignment");

  skip_grain_reassignment_steps =
    parameter_handler.get_integer("Time steps between grain reassignments");

  order_parameter_threshold =
    parameter_handler.get_double("Order parameter cutoff for grain identification");

  buffer_between_grains =
    parameter_handler.get_double("Buffer between grains before reassignment");
  if (buffer_between_grains < 0.0 && grain_remapping_activated)
    {
      std::cerr << "PRISMS-PF Error: If grain reassignment is activated, a "
                   "non-negative buffer distance must be given. See the 'Buffer "
                   "between grains before reassignment' entry in parameters.in.\n";
      abort();
    }

  const std::vector<std::string> variables_for_remapping_str =
    dealii::Utilities::split_string_list(
      parameter_handler.get("Order parameter fields for grain reassignment"));
  for (const auto &field : variables_for_remapping_str)
    {
      bool field_found = false;
      for (const auto &[index, variable] : variable_attributes.attributes)
        {
          if (boost::iequals(field, variable.name))
            {
              variables_for_remapping.push_back(index);
              field_found = true;
              break;
            }
        }
      if (!field_found && grain_remapping_activated)
        {
          std::cerr << "PRISMS-PF Error: Entries in the list of order "
                       "parameter fields used for grain reassignment must "
                       "match the variable names in equations.h.\n";
          std::cerr << field << "\n";
          abort();
        }
    }

  load_grain_structure = parameter_handler.get_bool("Load grain structure");
  load_vtk_file_type =
    parameter_handler.get("vtk file type"); // assign the vtk file type and getting it
                                            // ready to send to initialconditions.cc
  grain_structure_filename      = parameter_handler.get("Grain structure filename");
  grain_structure_variable_name = parameter_handler.get("Grain structure variable name");
  num_grain_smoothing_cycles    = parameter_handler.get_integer(
    "Number of smoothing cycles after grain structure loading");
  min_radius_for_loading_grains =
    parameter_handler.get_double("Minimum radius for loaded grains");
}

template <int dim>
void
userInputParameters<dim>::assign_boundary_condition_parameters(
  dealii::ParameterHandler &parameter_handler,
  variableAttributeLoader  &variable_attributes)
{
  // Load the boundary condition variables into list of BCs (where each element
  // of the vector is one component of one variable)
  std::vector<std::string> list_of_BCs;
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (variable.var_type == SCALAR)
        {
          std::string bc_text = "Boundary condition for variable ";
          bc_text.append(variable.name);
          list_of_BCs.push_back(parameter_handler.get(bc_text));
        }
      else
        {
          std::string bc_text = "Boundary condition for variable ";
          bc_text.append(variable.name);
          bc_text.append(", x component");
          list_of_BCs.push_back(parameter_handler.get(bc_text));

          bc_text = "Boundary condition for variable ";
          bc_text.append(variable.name);
          bc_text.append(", y component");
          list_of_BCs.push_back(parameter_handler.get(bc_text));

          if (dim > 2)
            {
              bc_text = "Boundary condition for variable ";
              bc_text.append(variable.name);
              bc_text.append(", z component");
              list_of_BCs.push_back(parameter_handler.get(bc_text));
            }
        }
    }

  /*----------------------
  |  Pinning point
  -----------------------*/
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      std::string pinning_text = "Pinning point: ";
      pinning_text.append(variable.name);
      parameter_handler.enter_subsection(pinning_text);

      // Skip if the default
      if (parameter_handler.get_double("x") == -1.0)
        {
          parameter_handler.leave_subsection();
          continue;
        }

      // Otherwise, fill out point
      if (dim == 2)
        {
          pinned_point[index] = dealii::Point<dim>(parameter_handler.get_double("x"),
                                                   parameter_handler.get_double("y"));
        }
      else
        {
          pinned_point[index] = dealii::Point<dim>(parameter_handler.get_double("x"),
                                                   parameter_handler.get_double("y"),
                                                   parameter_handler.get_double("z"));
        }
      parameter_handler.leave_subsection();
    }

  // Load the BC information from the strings into a varBCs object
  load_BC_list(list_of_BCs);
}

// NOLINTBEGIN(cppcoreguidelines-pro-type-member-init, hicpp-member-init)
template <int dim>
userInputParameters<dim>::userInputParameters(inputFileReader          &input_file_reader,
                                              dealii::ParameterHandler &parameter_handler,
                                              variableAttributeLoader variable_attributes)
{
  loadVariableAttributes(variable_attributes);

  // Spatial discretization
  assign_spatial_discretization_parameters(parameter_handler, variable_attributes);

  // Time stepping parameters
  assign_temporal_discretization_parameters(parameter_handler, variable_attributes);

  // Linear solver parameters
  assign_linear_solve_parameters(parameter_handler, variable_attributes);

  // Non-linear solver parameters
  assign_nonlinear_solve_parameters(parameter_handler, variable_attributes);

  // Output parameters
  assign_output_parameters(parameter_handler);

  // Initial condition parameters
  assign_load_initial_condition_parameters(parameter_handler);

  // Nucleation parameters
  assign_nucleation_parameters(parameter_handler, variable_attributes);

  // Grain remapping & vtk load-in parameters
  assign_grain_parameters(parameter_handler, variable_attributes);

  // Boundary conditions
  assign_boundary_condition_parameters(parameter_handler, variable_attributes);

  // Load the user-defined constants
  load_user_constants(input_file_reader, parameter_handler);
}

// NOLINTEND(cppcoreguidelines-pro-type-member-init, hicpp-member-init)

template <int dim>
void
userInputParameters<dim>::assign_boundary_conditions(
  std::vector<std::string> &boundary_condition_list,
  varBCs<dim>              &boundary_condition)
{
  // If there is only one boundary condition, copy it to have 2*dim copies.
  if (boundary_condition_list.size() == 1)
    {
      boundary_condition_list.resize(static_cast<size_t>(2 * dim),
                                     boundary_condition_list[0]);
    }

  // Assign the boundary condition into the varBCs<dim> object.
  for (unsigned int j = 0; j < (2 * dim); j++)
    {
      if (boost::iequals(boundary_condition_list[j], "NATURAL"))
        {
          boundary_condition.var_BC_type.push_back(NATURAL);
          boundary_condition.var_BC_val.push_back(0.0);
        }
      else if (boost::iequals(boundary_condition_list[j], "PERIODIC"))
        {
          boundary_condition.var_BC_type.push_back(PERIODIC);
          boundary_condition.var_BC_val.push_back(0.0);
        }
      else if (boost::iequals(boundary_condition_list[j], "NON_UNIFORM_DIRICHLET"))
        {
          boundary_condition.var_BC_type.push_back(NON_UNIFORM_DIRICHLET);
          boundary_condition.var_BC_val.push_back(0.0);
        }
      else if (boost::iequals(boundary_condition_list[j].substr(0, 9), "DIRICHLET"))
        {
          boundary_condition.var_BC_type.push_back(DIRICHLET);
          std::string dirichlet_val =
            boundary_condition_list[j].substr(10, boundary_condition_list[j].size());
          dirichlet_val = dealii::Utilities::trim(dirichlet_val);
          boundary_condition.var_BC_val.push_back(
            dealii::Utilities::string_to_double(dirichlet_val));
        }
      else if (boost::iequals(boundary_condition_list[j].substr(0, 7), "NEUMANN"))
        {
          boundary_condition.var_BC_type.push_back(NEUMANN);
          std::string neumann_val =
            boundary_condition_list[j].substr(8, boundary_condition_list[j].size());
          neumann_val = dealii::Utilities::trim(neumann_val);
          boundary_condition.var_BC_val.push_back(
            dealii::Utilities::string_to_double(neumann_val));
        }
      else
        {
          std::cout << boundary_condition_list[j].substr(0, 8) << "\n";
          std::cout << "Error: Boundary conditions specified improperly.\n";
          abort();
        }

      // If periodic BCs are used, ensure they are applied on both sides of
      // domain
      if (j % 2 == 0)
        {
          AssertThrow(boost::iequals(boundary_condition_list[j], "PERIODIC") ==
                        boost::iequals(boundary_condition_list[j + 1], "PERIODIC"),
                      dealii::ExcMessage(
                        std::string("Periodic boundary condition must be "
                                    "specified on both sides of domain")));
        }
    }
}

template <int dim>
void
userInputParameters<dim>::load_BC_list(const std::vector<std::string> &list_of_BCs)
{
  // Loop over the list of boundary conditions specified in parameters
  // and provided in the input list_of_BCs. Process the BCs and place
  // them into the vector BC_list
  std::vector<std::string> split_boundary_conditions;
  for (const auto &boundary_condition : list_of_BCs)
    {
      // Ensure all variables have BCs specified in parameters.prm
      AssertThrow(!boundary_condition.empty(),
                  dealii::ExcMessage(std::string("Boundary condition not specified.")));

      // Create object to store BCs
      varBCs<dim> newBC;

      // Split string of boundary conditions
      split_boundary_conditions =
        dealii::Utilities::split_string_list(boundary_condition);

      // Assign the boundaries into the object
      assign_boundary_conditions(split_boundary_conditions, newBC);

      // Append BCs for current field to total list
      BC_list.push_back(newBC);
    }
}

template <int dim>
std::vector<unsigned int>
userInputParameters<dim>::setTimeStepList(
  const std::string               &outputSpacingType,
  unsigned int                     numberOfOutputs,
  const std::vector<unsigned int> &userGivenTimeStepList)
{
  // Initialize timestep list
  std::vector<unsigned int> timeStepList;

  // The number of outputs cannot be greater than the number increments
  numberOfOutputs = std::min(numberOfOutputs, totalIncrements);

  // Prevent divide by zero in subsequent output types by returning the a vector where the
  // only entry is one greater than the number of increments. This way, we effectively
  // have no outputs. While this condition can be ignored for the LIST type, the user
  // should just ignore the parameter `set Number of outputs` and use the default value
  // of 10.
  if (numberOfOutputs == 0)
    {
      timeStepList.push_back(totalIncrements + 1);
      return timeStepList;
    }

  // Set output list for all the output list types
  if (outputSpacingType == "LIST")
    {
      timeStepList = userGivenTimeStepList;
    }
  else if (outputSpacingType == "EQUAL_SPACING")
    {
      for (unsigned int iter = 0; iter <= totalIncrements;
           iter += totalIncrements / numberOfOutputs)
        {
          timeStepList.push_back(iter);
        }
    }
  else if (outputSpacingType == "LOG_SPACING")
    {
      timeStepList.push_back(0);
      for (unsigned int output = 1; output <= numberOfOutputs; output++)
        {
          timeStepList.push_back(round(std::pow(static_cast<double>(totalIncrements),
                                                static_cast<double>(output) /
                                                  static_cast<double>(numberOfOutputs))));
        }
    }
  else if (outputSpacingType == "N_PER_DECADE")
    {
      AssertThrow(totalIncrements > 1,
                  dealii::ExcMessage(
                    std::string("PRISMS-PF Error: For n per decaded spaced outputs, "
                                "the number of increments must be greater than 1.")));

      timeStepList.push_back(0);
      timeStepList.push_back(1);
      for (unsigned int iter = 2; iter <= totalIncrements; iter++)
        {
          const unsigned int decade = std::ceil(std::log10(iter));
          const auto         step_size =
            static_cast<unsigned int>(std::pow(10, decade) / numberOfOutputs);
          if (iter % step_size == 0)
            {
              timeStepList.push_back(iter);
            }
        }
    }
  else
    {
      AssertThrow(false,
                  dealii::ExcMessage(
                    std::string("PRISMS-PF Error: Invalid output spacing type.")));
    }

  return timeStepList;
}

template <int dim>
void
userInputParameters<dim>::loadVariableAttributes(
  const variableAttributeLoader &variable_attributes)
{
  number_of_variables    = variable_attributes.attributes.size();
  pp_number_of_variables = variable_attributes.pp_attributes.size();
  // Load some nucleation parameters
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (variable.nucleating_variable)
        {
          nucleating_variable_indices.push_back(index);
        }
      if (variable.need_value_nucleation || variable.nucleating_variable)
        {
          nucleation_need_value.push_back(index);
        }
    }

  nucleation_occurs = !nucleating_variable_indices.empty();

  // Load variable information for calculating the RHS for explicit equations
  num_var_explicit_RHS = 0;
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (!static_cast<bool>(variable.eval_flags_explicit_RHS &
                             dealii::EvaluationFlags::nothing))
        {
          num_var_explicit_RHS++;
        }
    }
  varInfoListExplicitRHS.reserve(num_var_explicit_RHS);
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      variable_info varInfo {};

      varInfo.evaluation_flags = variable.eval_flags_explicit_RHS;

      varInfo.residual_flags = variable.eval_flags_residual_explicit_RHS;

      varInfo.global_var_index = index;

      varInfo.var_needed =
        !static_cast<bool>(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = variable.var_type == SCALAR;

      varInfoListExplicitRHS.push_back(varInfo);
    }

  // Load variable information for calculating the RHS for nonexplicit equations
  num_var_nonexplicit_RHS = 0;
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (!static_cast<bool>(variable.eval_flags_nonexplicit_RHS &
                             dealii::EvaluationFlags::nothing))
        {
          num_var_nonexplicit_RHS++;
        }
    }
  varInfoListNonexplicitRHS.reserve(num_var_nonexplicit_RHS);
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      variable_info varInfo {};

      varInfo.evaluation_flags = variable.eval_flags_nonexplicit_RHS;

      varInfo.residual_flags = variable.eval_flags_residual_nonexplicit_RHS;

      varInfo.global_var_index = index;

      varInfo.var_needed =
        !static_cast<bool>(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = variable.var_type == SCALAR;

      varInfoListNonexplicitRHS.push_back(varInfo);
    }

  // Load variable information for calculating the LHS
  num_var_LHS = 0;
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      if (!static_cast<bool>(variable.eval_flags_nonexplicit_LHS &
                             dealii::EvaluationFlags::nothing))
        {
          num_var_LHS++;
        }
    }

  varInfoListLHS.reserve(num_var_LHS);
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      variable_info varInfo {};

      varInfo.evaluation_flags = variable.eval_flags_nonexplicit_LHS;

      varInfo.residual_flags = variable.eval_flags_residual_nonexplicit_LHS;

      varInfo.global_var_index = index;

      varInfo.var_needed =
        !static_cast<bool>(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = variable.var_type == SCALAR;

      varInfoListLHS.push_back(varInfo);
    }

  varChangeInfoListLHS.reserve(num_var_LHS);
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      variable_info varInfo {};

      varInfo.evaluation_flags = variable.eval_flags_change_nonexplicit_LHS;

      // FOR NOW, TAKING THESE FROM THE VARIABLE ITSELF!!
      varInfo.residual_flags = variable.eval_flags_residual_nonexplicit_LHS;

      varInfo.global_var_index = index;

      varInfo.var_needed =
        !static_cast<bool>(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = variable.var_type == SCALAR;

      varChangeInfoListLHS.push_back(varInfo);
    }

  // Load variable information for postprocessing
  // First, the info list for the base field variables
  pp_baseVarInfoList.reserve(number_of_variables);
  for (const auto &[index, variable] : variable_attributes.attributes)
    {
      variable_info varInfo {};

      varInfo.evaluation_flags = variable.eval_flags_postprocess;

      varInfo.global_var_index = index;

      varInfo.var_needed =
        !static_cast<bool>(varInfo.evaluation_flags & dealii::EvaluationFlags::nothing);

      varInfo.is_scalar = variable.var_type == SCALAR;

      pp_baseVarInfoList.push_back(varInfo);
    }

  // Now load the information for the post-processing variables
  // Parameters for postprocessing

  postProcessingRequired = pp_number_of_variables > 0;

  num_integrated_fields = 0;
  for (const auto &[pp_index, pp_variable] : variable_attributes.pp_attributes)
    {
      if (pp_variable.calc_integral)
        {
          num_integrated_fields++;
          integrated_field_indices.push_back(pp_index);
        }
    }

  // The info list for the postprocessing field variables
  pp_varInfoList.reserve(pp_number_of_variables);
  for (const auto &[pp_index, pp_variable] : variable_attributes.pp_attributes)
    {
      variable_info varInfo {};

      varInfo.var_needed = true;

      varInfo.residual_flags = pp_variable.eval_flags_residual_postprocess;

      varInfo.global_var_index = pp_index;

      varInfo.is_scalar = pp_variable.var_type == SCALAR;

      pp_varInfoList.push_back(varInfo);
    }
}

template <int dim>
unsigned int
userInputParameters<dim>::compute_tensor_parentheses(
  const unsigned int              n_elements,
  const std::vector<std::string> &tensor_elements)
{
  unsigned int open_parentheses  = 0;
  unsigned int close_parentheses = 0;

  for (unsigned int element = 0; element < n_elements; element++)
    {
      for (const char c : tensor_elements.at(element))
        {
          if (c == '(')
            {
              ++open_parentheses;
            }
          else if (c == ')')
            {
              ++close_parentheses;
            }
        }
    }

  if (open_parentheses != close_parentheses)
    {
      std::cerr << "PRISMS-PF ERROR: User-defined elastic constant "
                   "list does not have the same number of open and "
                   "close parentheses.\n";
      abort();
    }

  return open_parentheses;
}

template <int dim>
void
userInputParameters<dim>::remove_parentheses(std::vector<std::string> &tensor_elements)
{
  for (std::string &element : tensor_elements)
    {
      element.erase(std::remove(element.begin(), element.end(), '('), element.end());
      element.erase(std::remove(element.begin(), element.end(), ')'), element.end());
    }
}

template <int dim>
dealii::Tensor<1, dim>
userInputParameters<dim>::compute_rank_1_tensor_constant(
  const unsigned int       n_elements,
  std::vector<std::string> tensor_elements)
{
  AssertThrow(n_elements > 1 && n_elements < 4,
              dealii::ExcMessage("PRISMS-PF Error: The columns in user-defined constant "
                                 "tensors cannot be longer than 3 elements (internally "
                                 "truncated to the number of dimensions)."));

  dealii::Tensor<1, dim> temp;
  for (unsigned int i = 0; i < dim; i++)
    {
      temp[i] = dealii::Utilities::string_to_double(tensor_elements.at(i));
    }

  return temp;
}

template <int dim>
dealii::Tensor<2, dim>
userInputParameters<dim>::compute_rank_2_tensor_constant(
  const unsigned int       n_elements,
  std::vector<std::string> tensor_elements)
{
  unsigned int row_length = 0;
  if (n_elements == 4)
    {
      AssertThrow(dim < 3,
                  dealii::ExcMessage(
                    "PRISMS-PF ERROR: User-defined constant tensor does not have "
                    "enough elements. For 3D calculations matrices must be 3x3."));

      row_length = 2;
    }
  else if (n_elements == 9)
    {
      row_length = 3;
    }
  else
    {
      AssertThrow(false,
                  dealii::ExcMessage("PRISMS-PF ERROR: User-defined constant tensor does "
                                     "not have the correct number of elements, matrices "
                                     "must be 2x2 or 3x3."));
    }

  dealii::Tensor<2, dim> temp;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          temp[i][j] =
            dealii::Utilities::string_to_double(tensor_elements.at(i * row_length + j));
        }
    }

  return temp;
}

template <int dim>
void
userInputParameters<dim>::assign_user_constant(
  std::vector<std::string> &model_constants_strings)
{
  // Ensure that the input includes a value and a type
  AssertThrow(model_constants_strings.size() > 1,
              dealii::ExcMessage("PRISMS-PF Error: At least two fields are required for "
                                 "user-defined variables (value and type)."));

  std::vector<std::string> model_constants_type_strings =
    dealii::Utilities::split_string_list(model_constants_strings.at(
                                           model_constants_strings.size() - 1),
                                         ' ');

  if (model_constants_strings.size() == 2)
    {
      assign_primitive_user_constant(model_constants_strings);
    }
  else
    {
      if (boost::iequals(model_constants_type_strings.at(0), "tensor"))
        {
          const unsigned int n_elements = model_constants_strings.size() - 1;

          const unsigned int open_parentheses =
            compute_tensor_parentheses(n_elements, model_constants_strings);
          remove_parentheses(model_constants_strings);

          // Rank 1 tensor
          if (open_parentheses < 3)
            {
              model_constants.push_back(
                compute_rank_1_tensor_constant(n_elements, model_constants_strings));
            }
          // Rank 2 tensor
          else if (open_parentheses < 5)
            {
              model_constants.push_back(
                compute_rank_2_tensor_constant(n_elements, model_constants_strings));
            }
        }
      else if (boost::iequals(model_constants_type_strings.at(1), "elastic") &&
               boost::iequals(model_constants_type_strings.at(2), "constants"))
        {
          const unsigned int n_elements = model_constants_strings.size() - 1;

          remove_parentheses(model_constants_strings);

          // Load in the elastic constants as a vector
          std::vector<double> temp_elastic_constants;
          for (unsigned int i = 0; i < n_elements; i++)
            {
              temp_elastic_constants.push_back(
                dealii::Utilities::string_to_double(model_constants_strings.at(i)));
            }

          const std::string elastic_const_symmetry = model_constants_type_strings.at(0);
          dealii::Tensor<2, 2 *dim - 1 + dim / 3> temp =
            get_Cij_tensor(temp_elastic_constants, elastic_const_symmetry);
          model_constants.push_back(temp);
        }
      else
        {
          AssertThrow(false,
                      dealii::ExcMessage(
                        "PRISMS-PF ERROR: Only user-defined constant tensors may "
                        "have multiple elements."));
        }
    }
}

template <int dim>
void
userInputParameters<dim>::assign_primitive_user_constant(
  std::vector<std::string> &model_constants_strings)
{
  std::vector<std::string> model_constants_type_strings =
    dealii::Utilities::split_string_list(model_constants_strings.at(
                                           model_constants_strings.size() - 1),
                                         ' ');

  if (boost::iequals(model_constants_type_strings.at(0), "double"))
    {
      model_constants.push_back(
        dealii::Utilities::string_to_double(model_constants_strings.at(0)));
    }
  else if (boost::iequals(model_constants_type_strings.at(0), "int"))
    {
      model_constants.push_back(
        dealii::Utilities::string_to_int(model_constants_strings.at(0)));
    }
  else if (boost::iequals(model_constants_type_strings.at(0), "bool"))
    {
      bool temp = boost::iequals(model_constants_strings.at(0), "true");
      model_constants.push_back(temp);
    }
  else
    {
      AssertThrow(false,
                  dealii::ExcMessage(
                    "PRISMS-PF Error: The type for user-defined variables must be "
                    "`double`, `int`, `bool`, `tensor`, or `elastic constants`."));
    }
}

template <int dim>
void
userInputParameters<dim>::load_user_constants(inputFileReader          &input_file_reader,
                                              dealii::ParameterHandler &parameter_handler)
{
  const unsigned int number_of_constants = input_file_reader.num_constants;

  for (unsigned int i = 0; i < input_file_reader.model_constant_names.size(); i++)
    {
      model_constant_name_map[input_file_reader.model_constant_names[i]] = i;
    }

  for (unsigned int i = 0; i < number_of_constants; i++)
    {
      std::string constants_text = "Model constant ";
      constants_text.append(input_file_reader.model_constant_names[i]);

      std::vector<std::string> model_constants_strings =
        dealii::Utilities::split_string_list(parameter_handler.get(constants_text));

      assign_user_constant(model_constants_strings);
    }
}

template <int dim>
dealii::Tensor<2, 2 * dim - 1 + dim / 3>
userInputParameters<dim>::get_Cij_tensor(std::vector<double> elastic_constants,
                                         const std::string  &elastic_const_symmetry) const
{
  // First set the material model
  elasticityModel mat_model = ISOTROPIC;
  if (elastic_const_symmetry == "isotropic")
    {
      mat_model = ISOTROPIC;
    }
  else if (elastic_const_symmetry == "transverse")
    {
      mat_model = TRANSVERSE;
    }
  else if (elastic_const_symmetry == "orthotropic")
    {
      mat_model = ORTHOTROPIC;
    }
  else if (elastic_const_symmetry == "anisotropic")
    {
      mat_model = ANISOTROPIC;
    }
  else
    {
      // Should change to an exception
      std::cerr << "Elastic material model is invalid, please use isotropic, "
                   "transverse, orthotropic, or anisotropic\n";
    }

  // If the material model is anisotropic for a 2D calculation but the elastic
  // constants are given for a 3D calculation, change the elastic constant
  // vector to the 2D form
  if ((mat_model == ANISOTROPIC) && (dim == 2) && elastic_constants.size() == 21)
    {
      std::vector<double> elastic_constants_temp = elastic_constants;
      elastic_constants.clear();
      const std::vector<unsigned int> indices_2D = {0, 1, 5, 6, 10, 14};
      for (const auto &index : indices_2D)
        {
          elastic_constants.push_back(elastic_constants_temp.at(index));
        }
    }

  dealii::ConditionalOStream pcout(std::cout,
                                   dealii::Utilities::MPI::this_mpi_process(
                                     MPI_COMM_WORLD) == 0);

  return getCIJMatrix(mat_model, elastic_constants, pcout);
}

template <int dim>
dealii::Tensor<2, 2 * dim - 1 + dim / 3>
userInputParameters<dim>::getCIJMatrix(const elasticityModel       model,
                                       const std::vector<double>  &constants,
                                       dealii::ConditionalOStream &pcout) const
{
  // CIJ.fill(0.0);
  dealii::Tensor<2, 2 * dim - 1 + dim / 3> CIJ;

  pcout << "Reading material model:";
  switch (dim)
    {
      case 1:
        {
          pcout << " 1D ";
          // 1D models
          switch (model)
            {
              case ISOTROPIC:
                {
                  pcout << " ISOTROPIC \n";
                  CIJ[0][0] = constants[0];
                  break;
                }
              default:
                {
                  std::cout << "\nelasticityModels: Supported models in 1D - "
                               "ISOTROPIC\n";
                  std::cout << "See /src/elasticityModels.h\n";
                  exit(-1);
                }
            }
          break;
        }
      case 2:
        {
          pcout << " 2D ";
          // 2D models
          switch (model)
            {
              case ISOTROPIC:
                {
                  pcout << " ISOTROPIC \n";
                  const double E      = constants[0];
                  const double nu     = constants[1];
                  const double mu     = E / (2 * (1 + nu));
                  const double lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
                  CIJ[0][0]           = lambda + 2 * mu;
                  CIJ[1][1]           = lambda + 2 * mu;
                  CIJ[2][2]           = mu;
                  CIJ[0][1] = CIJ[1][0] = lambda;
                  break;
                }
              case ANISOTROPIC:
                {
                  pcout << " ANISOTROPIC \n";
                  CIJ[0][0] = constants[0];             // C11
                  CIJ[1][1] = constants[1];             // C22
                  CIJ[2][2] = constants[2];             // C33
                  CIJ[0][1] = CIJ[1][0] = constants[3]; // C12
                  CIJ[0][2] = CIJ[2][0] = constants[4]; // C13
                  CIJ[1][2] = CIJ[2][1] = constants[5]; // C23
                  break;
                }
              default:
                {
                  std::cout << "\nelasticityModels: Supported models in 2D - "
                               "ISOTROPIC/ANISOTROPIC\n";
                  std::cout << "See /src/elasticityModels.h\n";
                  exit(-1);
                }
            }
          break;
        }
      case 3:
        {
          pcout << " 3D ";
          // 3D models
          switch (model)
            {
              case ISOTROPIC:
                {
                  pcout << " ISOTROPIC \n";
                  const double E      = constants[0];
                  const double nu     = constants[1];
                  const double mu     = E / (2 * (1 + nu));
                  const double lambda = nu * E / ((1 + nu) * (1 - 2 * nu));
                  CIJ[0][0]           = lambda + 2 * mu;
                  CIJ[1][1]           = lambda + 2 * mu;
                  CIJ[2][2]           = lambda + 2 * mu;
                  CIJ[3][3]           = mu;
                  CIJ[4][4]           = mu;
                  CIJ[5][5]           = mu;
                  CIJ[0][1] = CIJ[1][0] = lambda;
                  CIJ[0][2] = CIJ[2][0] = lambda;
                  CIJ[1][2] = CIJ[2][1] = lambda;
                  break;
                }
              case TRANSVERSE:
                {
                  pcout << " TRANSVERSE \n";
                  CIJ[0][0] = constants[0];                        // C11
                  CIJ[1][1] = constants[0];                        // C11
                  CIJ[2][2] = constants[1];                        // C33
                  CIJ[3][3] = constants[2];                        // C44
                  CIJ[4][4] = constants[2];                        // C44
                  CIJ[5][5] = (constants[0] - constants[3]) / 2.0; //(C11-C12)/2
                  CIJ[0][1] = CIJ[1][0] = constants[3];            // C12
                  CIJ[0][2] = CIJ[2][0] = constants[4];            // C13
                  CIJ[1][2] = CIJ[2][1] = constants[4];            // C13
                  break;
                }
              case ORTHOTROPIC:
                {
                  pcout << " ORTHOTROPIC \n";
                  CIJ[0][0] = constants[0];             // C11
                  CIJ[1][1] = constants[1];             // C22
                  CIJ[2][2] = constants[2];             // C33
                  CIJ[3][3] = constants[3];             // C44
                  CIJ[4][4] = constants[4];             // C55
                  CIJ[5][5] = constants[5];             // C66
                  CIJ[0][1] = CIJ[1][0] = constants[6]; // C12
                  CIJ[0][2] = CIJ[2][0] = constants[7]; // C13
                  CIJ[1][2] = CIJ[2][1] = constants[8]; // C23
                  break;
                }
              case ANISOTROPIC:
                {
                  pcout << " ANISOTROPIC \n";
                  CIJ[0][0] = constants[0];              // C11
                  CIJ[1][1] = constants[1];              // C22
                  CIJ[2][2] = constants[2];              // C33
                  CIJ[3][3] = constants[3];              // C44
                  CIJ[4][4] = constants[4];              // C55
                  CIJ[5][5] = constants[5];              // C66
                  CIJ[0][1] = CIJ[1][0] = constants[6];  // C12
                  CIJ[0][2] = CIJ[2][0] = constants[7];  // C13
                  CIJ[0][3] = CIJ[3][0] = constants[8];  // C14
                  CIJ[0][4] = CIJ[4][0] = constants[9];  // C15
                  CIJ[0][5] = CIJ[5][0] = constants[10]; // C16
                  CIJ[1][2] = CIJ[2][1] = constants[11]; // C23
                  CIJ[1][3] = CIJ[3][1] = constants[12]; // C24
                  CIJ[1][4] = CIJ[4][1] = constants[13]; // C25
                  CIJ[1][5] = CIJ[5][1] = constants[14]; // C26
                  CIJ[2][3] = CIJ[3][2] = constants[15]; // C34
                  CIJ[2][4] = CIJ[4][2] = constants[16]; // C35
                  CIJ[2][5] = CIJ[5][2] = constants[17]; // C36
                  CIJ[3][4] = CIJ[4][3] = constants[18]; // C45
                  CIJ[3][5] = CIJ[5][3] = constants[19]; // C46
                  CIJ[4][5] = CIJ[5][4] = constants[20]; // C56
                  break;
                }
              default:
                {
                  std::cout << "\nelasticityModels: Supported models in 3D - "
                               "ISOTROPIC/TRANSVERSE/ORTHOTROPIC/ANISOTROPIC\n";
                  std::cout << "See /src/elasticityModels.h\n";
                  exit(-1);
                }
            }
          break;
        }
      default:
        {
          std::cout << "\nelasticityModels: DIM is not 1/2/3\n";
          exit(-1);
        }
    }
  // print CIJ to terminal
  pcout << "Elasticity matrix (Voigt notation):\n";
  constexpr unsigned int voight_matrix_size = 2 * dim - 1 + dim / 3;
  for (unsigned int i = 0; i < voight_matrix_size; i++)
    {
      for (unsigned int j = 0; j < voight_matrix_size; j++)
        {
          pcout << std::setw(8) << std::setprecision(3) << std::scientific << CIJ[i][j]
                << " ";
        }
      pcout << "\n";
    }
  pcout << "\n";
  return CIJ;
}

template class userInputParameters<2>;
template class userInputParameters<3>;
