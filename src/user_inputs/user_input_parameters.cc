// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <boost/algorithm/string/predicate.hpp>

#include <prismspf/core/exceptions.h>
#include <prismspf/core/refinement_criterion.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/input_file_reader.h>
#include <prismspf/user_inputs/linear_solve_parameters.h>
#include <prismspf/user_inputs/nonlinear_solve_parameters.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <climits>
#include <cmath>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
UserInputParameters<dim>::UserInputParameters(InputFileReader          &input_file_reader,
                                              dealii::ParameterHandler &parameter_handler)
  : var_attributes(input_file_reader.get_var_attributes())
{
  // Assign the parameters to the appropriate data structures
  assign_spatial_discretization_parameters(parameter_handler);
  assign_temporal_discretization_parameters(parameter_handler);
  assign_linear_solve_parameters(parameter_handler);
  assign_nonlinear_solve_parameters(parameter_handler);
  assign_output_parameters(parameter_handler);
  assign_checkpoint_parameters(parameter_handler);
  assign_boundary_parameters(parameter_handler);
  assign_load_initial_condition_parameters(parameter_handler);
  load_model_constants(input_file_reader, parameter_handler);

  // Perform and postprocessing of user inputs and run checks
  spatial_discretization.postprocess_and_validate();
  temporal_discretization.postprocess_and_validate(var_attributes);
  linear_solve_parameters.postprocess_and_validate();
  nonlinear_solve_parameters.postprocess_and_validate();
  output_parameters.postprocess_and_validate(temporal_discretization);
  checkpoint_parameters.postprocess_and_validate(temporal_discretization);
  boundary_parameters.postprocess_and_validate(var_attributes);
  load_ic_parameters.postprocess_and_validate();

  // Print all the parameters to summary.log
  spatial_discretization.print_parameter_summary();
  temporal_discretization.print_parameter_summary();
  linear_solve_parameters.print_parameter_summary();
  nonlinear_solve_parameters.print_parameter_summary();
  output_parameters.print_parameter_summary();
  checkpoint_parameters.print_parameter_summary();
  boundary_parameters.print_parameter_summary();
  load_ic_parameters.print_parameter_summary();
  user_constants.print();
}

template <unsigned int dim>
void
UserInputParameters<dim>::assign_spatial_discretization_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  parameter_handler.enter_subsection("Rectangular mesh");
  {
    std::vector<std::string> axis_labels = {"x", "y", "z"};
    for (unsigned int i = 0; i < dim; ++i)
      {
        spatial_discretization.set_size(i,
                                        parameter_handler.get_double(axis_labels[i] +
                                                                     " size"));
        spatial_discretization.set_subdivisions(i,
                                                parameter_handler.get_integer(
                                                  axis_labels[i] + " subdivisions"));
      }
  }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Spherical mesh");
  {
    spatial_discretization.set_radius(parameter_handler.get_double("radius"));
  }
  parameter_handler.leave_subsection();

  spatial_discretization.set_global_refinement(
    parameter_handler.get_integer("global refinement"));

  spatial_discretization.set_degree(parameter_handler.get_integer("degree"));

  spatial_discretization.set_has_adaptivity(
    parameter_handler.get_bool("mesh adaptivity"));

  spatial_discretization.set_remeshing_period(
    parameter_handler.get_integer("remeshing period"));

  spatial_discretization.set_max_refinement(
    parameter_handler.get_integer("max refinement"));
  spatial_discretization.set_min_refinement(
    parameter_handler.get_integer("min refinement"));

  for (const auto &[index, variable] : var_attributes)
    {
      std::string subsection_text = "refinement criterion: ";
      subsection_text.append(variable.get_name());
      parameter_handler.enter_subsection(subsection_text);
      {
        const std::string crit_type_string = parameter_handler.get("type");
        if (!boost::iequals(crit_type_string, "none"))
          {
            GridRefinement::RefinementCriterion new_criterion(
              GridRefinement::RefinementFlags::Nothing,
              parameter_handler.get_double("value lower bound"),
              parameter_handler.get_double("value upper bound"),
              parameter_handler.get_double("gradient magnitude lower bound"));

            if (boost::iequals(crit_type_string, "value"))
              {
                new_criterion.set_criterion(GridRefinement::RefinementFlags::Value);
              }
            else if (boost::iequals(crit_type_string, "gradient"))
              {
                new_criterion.set_criterion(GridRefinement::RefinementFlags::Gradient);
              }
            else if (boost::iequals(crit_type_string, "value_and_gradient"))
              {
                new_criterion.set_criterion(GridRefinement::RefinementFlags::Value |
                                            GridRefinement::RefinementFlags::Gradient);
              }
            else
              {
                AssertThrow(false, UnreachableCode());
              }
            spatial_discretization.add_refinement_criteria(new_criterion);
          }
      }
      parameter_handler.leave_subsection();
    }
}

template <unsigned int dim>
void
UserInputParameters<dim>::assign_temporal_discretization_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  temporal_discretization.set_timestep(parameter_handler.get_double("time step"));
  temporal_discretization.set_final_time(parameter_handler.get_double("end time"));
  temporal_discretization.set_total_increments(
    static_cast<unsigned int>(parameter_handler.get_integer("number steps")));
}

template <unsigned int dim>
void
UserInputParameters<dim>::assign_output_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  parameter_handler.enter_subsection("output");
  {
    output_parameters.set_file_name(parameter_handler.get("file name"));
    output_parameters.set_file_type(parameter_handler.get("file type"));
    output_parameters.set_patch_subdivisions(
      parameter_handler.get_integer("subdivisions"));
    output_parameters.set_output_condition(parameter_handler.get("condition"));
    output_parameters.set_user_output_list(dealii::Utilities::string_to_int(
      dealii::Utilities::split_string_list(parameter_handler.get("list"))));
    output_parameters.set_n_outputs(
      static_cast<unsigned int>(parameter_handler.get_integer("number")));
    output_parameters.set_print_output_period(
      parameter_handler.get_integer("print step period"));
    output_parameters.set_print_timing_with_output(
      parameter_handler.get_bool("timing information with output"));
  }
  parameter_handler.leave_subsection();
}

template <unsigned int dim>
void
UserInputParameters<dim>::assign_checkpoint_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  parameter_handler.enter_subsection("checkpoints");
  {
    checkpoint_parameters.set_load_from_checkpoint(
      parameter_handler.get_bool("load from checkpoint"));
    checkpoint_parameters.set_condition(parameter_handler.get("condition"));
    checkpoint_parameters.set_user_checkpoint_list(dealii::Utilities::string_to_int(
      dealii::Utilities::split_string_list(parameter_handler.get("list"))));
    checkpoint_parameters.set_n_checkpoints(
      static_cast<unsigned int>(parameter_handler.get_integer("number")));
  }
  parameter_handler.leave_subsection();
}

template <unsigned int dim>
void
UserInputParameters<dim>::assign_boundary_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  std::vector<std::string> axis_labels = {"x", "y", "z"};

  // Assign the normal boundary parameters
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.is_postprocess())
        {
          continue;
        }

      const unsigned int n_components =
        (variable.get_field_type() == FieldType::Scalar) ? 1 : dim;
      for (unsigned int i = 0; i < n_components; i++)
        {
          std::string bc_text = "boundary condition for " + variable.get_name();
          if (variable.get_field_type() != FieldType::Scalar)
            {
              bc_text += ", " + axis_labels[i] + " component";
            }
          boundary_parameters
            .set_boundary_condition_string(parameter_handler.get(bc_text), index, i);
        }
    }

  // Assign any pinning points
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.is_postprocess())
        {
          continue;
        }
      std::string pinning_text = "pinning point for ";
      pinning_text.append(variable.get_name());
      parameter_handler.enter_subsection(pinning_text);

      const std::string value_key =
        (variable.get_field_type() == FieldType::Scalar) ? "value" : "x value";

      // Skip if the value is the default INT_MAX
      if (parameter_handler.get_double(value_key) == INT_MAX)
        {
          parameter_handler.leave_subsection();
          continue;
        }

      // Fill out the point and value
      dealii::Point<dim> point;
      for (unsigned int i = 0; i < dim; ++i)
        {
          point[i] = parameter_handler.get_double(axis_labels[i]);
        }

      if (variable.get_field_type() == FieldType::Scalar)
        {
          boundary_parameters.set_pinned_point(parameter_handler.get_double("value"),
                                               point,
                                               index);
        }
      else
        {
          std::vector<double> value(dim);
          for (unsigned int i = 0; i < dim; ++i)
            {
              value[i] = parameter_handler.get_double(axis_labels[i] + " value");
            }
          boundary_parameters.set_pinned_point(value, point, index);
        }
      parameter_handler.leave_subsection();
    }
}

template <unsigned int dim>
void
UserInputParameters<dim>::assign_linear_solve_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.get_pde_type() == PDEType::TimeIndependent ||
          variable.get_pde_type() == PDEType::ImplicitTimeDependent)
        {
          std::string subsection_text = "linear solver parameters: ";
          subsection_text.append(variable.get_name());
          parameter_handler.enter_subsection(subsection_text);

          LinearSolverParameters linear_solver_parameters;

          // Set the tolerance type
          const std::string type_string = parameter_handler.get("tolerance type");
          if (boost::iequals(type_string, "AbsoluteResidual"))
            {
              linear_solver_parameters.tolerance_type =
                SolverToleranceType::AbsoluteResidual;
            }
          else if (boost::iequals(type_string, "RelativeResidualChange"))
            {
              linear_solver_parameters.tolerance_type =
                SolverToleranceType::RelativeResidualChange;
            }
          else
            {
              AssertThrow(false, UnreachableCode());
            }

          // Set the tolerance value
          linear_solver_parameters.tolerance =
            parameter_handler.get_double("tolerance value");

          // Set the maximum number of iterations
          linear_solver_parameters.max_iterations =
            parameter_handler.get_integer("max iterations");

          // Set preconditioner type and related parameters
          linear_solver_parameters.preconditioner =
            boost::iequals(parameter_handler.get("preconditioner type"), "GMG")
              ? PreconditionerType::GMG
              : PreconditionerType::None;

          linear_solver_parameters.smoothing_range =
            parameter_handler.get_double("smoothing range");

          linear_solver_parameters.smoother_degree =
            parameter_handler.get_integer("smoother degree");

          linear_solver_parameters.eig_cg_n_iterations =
            parameter_handler.get_integer("eigenvalue cg iterations");

          linear_solver_parameters.min_mg_level =
            parameter_handler.get_integer("min mg level");

          linear_solve_parameters.set_linear_solve_parameters(index,
                                                              linear_solver_parameters);

          parameter_handler.leave_subsection();
        }
    }
}

template <unsigned int dim>
void
UserInputParameters<dim>::assign_nonlinear_solve_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.get_field_solve_type() == FieldSolveType::NonexplicitSelfnonlinear ||
          variable.get_field_solve_type() == FieldSolveType::NonexplicitCononlinear)
        {
          std::string subsection_text = "nonlinear solver parameters: ";
          subsection_text.append(variable.get_name());
          parameter_handler.enter_subsection(subsection_text);

          NonlinearSolverParameters nonlinear_solver_parameters;
          nonlinear_solver_parameters.max_iterations =
            parameter_handler.get_integer("max iterations");
          nonlinear_solver_parameters.step_length =
            parameter_handler.get_double("step size");
          nonlinear_solve_parameters
            .set_nonlinear_solve_parameters(index, nonlinear_solver_parameters);

          // TODO (landinjm): Implement backtracking line search

          parameter_handler.leave_subsection();
        }
    }
}

template <unsigned int dim>
void
UserInputParameters<dim>::assign_load_initial_condition_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  load_ic_parameters.set_read_initial_conditions_from_file(
    parameter_handler.get_bool("read initial conditions from file"));

  for (unsigned int i = 0; i < 8; i++)
    {
      parameter_handler.enter_subsection("initial condition file " + std::to_string(i));
      {
        // Check if the file is specified
        if (parameter_handler.get("file name") != "")
          {
            // Create the LoadICFile object
            InitialConditionFile ic_file;
            ic_file.filename            = parameter_handler.get("file name");
            ic_file.file_extension      = parameter_handler.get("file extension");
            ic_file.grid_type           = parameter_handler.get("grid type");
            ic_file.file_variable_names = dealii::Utilities::split_string_list(
              parameter_handler.get("file variable names"));
            ic_file.simulation_variable_names = dealii::Utilities::split_string_list(
              parameter_handler.get("simulation variable names"));

            load_ic_parameters.add_ic_file(ic_file);
          }
      }
      parameter_handler.leave_subsection();
    }
}

template <unsigned int dim>
void
UserInputParameters<dim>::load_model_constants(
  const InputFileReader    &input_file_reader,
  dealii::ParameterHandler &parameter_handler)
{
  for (const std::string &constant_name : input_file_reader.get_model_constant_names())
    {
      std::string constants_text = "Model constant ";
      constants_text.append(constant_name);

      std::vector<std::string> model_constants_strings =
        dealii::Utilities::split_string_list(parameter_handler.get(constants_text));

      user_constants.add_user_constant(constant_name, model_constants_strings);
    }
}

INSTANTIATE_UNI_TEMPLATE(UserInputParameters)

PRISMS_PF_END_NAMESPACE
