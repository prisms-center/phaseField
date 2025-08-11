#include <deal.II/base/exceptions.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/utilities.h>

#include <boost/algorithm/string/predicate.hpp>

#include <prismspf/config.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/refinement_criterion.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/user_inputs/input_file_reader.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
userInputParameters<dim>::userInputParameters(inputFileReader          &input_file_reader,
                                              dealii::ParameterHandler &parameter_handler)
  : var_attributes(input_file_reader.var_attributes)
{
  // Assign the parameters to the appropriate data structures
  assign_spatial_discretization_parameters(parameter_handler);
  assign_temporal_discretization_parameters(parameter_handler);
  assign_linear_solve_parameters(parameter_handler);
  assign_nonlinear_solve_parameters(parameter_handler);
  assign_output_parameters(parameter_handler);
  assign_checkpoint_parameters(parameter_handler);
  assign_boundary_parameters(parameter_handler);
  load_model_constants(input_file_reader, parameter_handler);

  // Perform and postprocessing of user inputs and run checks
  spatial_discretization.postprocess_and_validate();
  temporal_discretization.postprocess_and_validate(var_attributes);
  linear_solve_parameters.postprocess_and_validate();
  nonlinear_solve_parameters.postprocess_and_validate();
  output_parameters.postprocess_and_validate(temporal_discretization);
  checkpoint_parameters.postprocess_and_validate(temporal_discretization);
  boundary_parameters.postprocess_and_validate(var_attributes);

  // Print all the parameters to summary.log
  spatial_discretization.print_parameter_summary();
  temporal_discretization.print_parameter_summary();
  linear_solve_parameters.print_parameter_summary();
  nonlinear_solve_parameters.print_parameter_summary();
  output_parameters.print_parameter_summary();
  checkpoint_parameters.print_parameter_summary();
  boundary_parameters.print_parameter_summary();
}

template <int dim>
void
userInputParameters<dim>::assign_spatial_discretization_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  parameter_handler.enter_subsection("rectangular mesh");
  {
    std::vector<std::string> axis_labels = {"x", "y", "z"};
    for (unsigned int i = 0; i < dim; ++i)
      {
        spatial_discretization.size[i] =
          parameter_handler.get_double(axis_labels[i] + " size");
        spatial_discretization.subdivisions[i] =
          parameter_handler.get_integer(axis_labels[i] + " subdivisions");
      }
  }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("spherical mesh");
  {
    spatial_discretization.radius = parameter_handler.get_double("radius");
  }
  parameter_handler.leave_subsection();

  spatial_discretization.global_refinement =
    parameter_handler.get_integer("global refinement");

  spatial_discretization.degree = parameter_handler.get_integer("degree");

  spatial_discretization.has_adaptivity = parameter_handler.get_bool("mesh adaptivity");

  spatial_discretization.remeshing_period =
    parameter_handler.get_integer("remeshing period");

  spatial_discretization.max_refinement = parameter_handler.get_integer("max refinement");
  spatial_discretization.min_refinement = parameter_handler.get_integer("min refinement");

  for (const auto &[index, variable] : var_attributes)
    {
      std::string subsection_text = "refinement criterion: ";
      subsection_text.append(variable.name);
      parameter_handler.enter_subsection(subsection_text);
      {
        const std::string crit_type_string = parameter_handler.get("type");
        if (!boost::iequals(crit_type_string, "none"))
          {
            RefinementCriterion new_criterion;
            new_criterion.variable_index = index;
            new_criterion.variable_name  = variable.name;
            if (boost::iequals(crit_type_string, "value"))
              {
                new_criterion.criterion_type = criterion_value;
                new_criterion.value_lower_bound =
                  parameter_handler.get_double("value lower bound");
                new_criterion.value_upper_bound =
                  parameter_handler.get_double("value upper bound");
              }
            else if (boost::iequals(crit_type_string, "gradient"))
              {
                new_criterion.criterion_type = criterion_gradient;
                new_criterion.gradient_lower_bound =
                  parameter_handler.get_double("gradient magnitude lower bound");
              }
            else if (boost::iequals(crit_type_string, "value_and_gradient"))
              {
                new_criterion.criterion_type = criterion_value | criterion_gradient;
                new_criterion.value_lower_bound =
                  parameter_handler.get_double("value lower bound");
                new_criterion.value_upper_bound =
                  parameter_handler.get_double("value upper bound");
                new_criterion.gradient_lower_bound =
                  parameter_handler.get_double("gradient magnitude lower bound");
              }
            else
              {
                AssertThrow(false, UnreachableCode());
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
  temporal_discretization.dt         = parameter_handler.get_double("time step");
  temporal_discretization.final_time = parameter_handler.get_double("end time");
  temporal_discretization.total_increments =
    static_cast<unsigned int>(parameter_handler.get_integer("number steps"));
}

template <int dim>
void
userInputParameters<dim>::assign_output_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  parameter_handler.enter_subsection("output");
  {
    output_parameters.file_name = parameter_handler.get("file name");
    output_parameters.file_type = parameter_handler.get("file type");
    output_parameters.output_per_process =
      parameter_handler.get_bool("separate files per process");
    output_parameters.condition        = parameter_handler.get("condition");
    output_parameters.user_output_list = dealii::Utilities::string_to_int(
      dealii::Utilities::split_string_list(parameter_handler.get("list")));
    output_parameters.n_outputs =
      static_cast<unsigned int>(parameter_handler.get_integer("number"));
    output_parameters.print_output_period =
      parameter_handler.get_integer("print step period");
    output_parameters.print_timing_with_output =
      parameter_handler.get_bool("timing information with output");
  }
  parameter_handler.leave_subsection();
}

template <int dim>
void
userInputParameters<dim>::assign_checkpoint_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  parameter_handler.enter_subsection("checkpoints");
  {
    checkpoint_parameters.load_from_checkpoint =
      parameter_handler.get_bool("load from checkpoint");
    checkpoint_parameters.condition            = parameter_handler.get("condition");
    checkpoint_parameters.user_checkpoint_list = dealii::Utilities::string_to_int(
      dealii::Utilities::split_string_list(parameter_handler.get("list")));
    checkpoint_parameters.n_checkpoints =
      static_cast<unsigned int>(parameter_handler.get_integer("number"));
  }
  parameter_handler.leave_subsection();
}

template <int dim>
void
userInputParameters<dim>::assign_boundary_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.is_postprocess)
        {
          continue;
        }
      if (variable.field_type == SCALAR)
        {
          std::string bc_text = "boundary condition for ";
          bc_text.append(variable.name);
          boundary_parameters.BC_list[index].emplace(0, parameter_handler.get(bc_text));
        }
      else
        {
          std::vector<std::string> axis_labels = {"x", "y", "z"};
          for (unsigned int i = 0; i < dim; i++)
            {
              std::string bc_text = "boundary condition for ";
              bc_text.append(variable.name + ", " + axis_labels[i] + " component");
              boundary_parameters.BC_list[index].emplace(i,
                                                         parameter_handler.get(bc_text));
            }
        }
    }

  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.is_postprocess)
        {
          continue;
        }
      std::string pinning_text = "pinning point for ";
      pinning_text.append(variable.name);
      parameter_handler.enter_subsection(pinning_text);
      if (variable.field_type == SCALAR)
        {
          // Skip if the value is the default INT_MAX
          if (parameter_handler.get_double("value") == 2147483647)
            {
              parameter_handler.leave_subsection();
              continue;
            }
          // Otherwise, fill out point and value
          std::vector<std::string> axis_labels = {"x", "y", "z"};
          dealii::Point<dim>       point;
          for (unsigned int i = 0; i < dim; ++i)
            {
              point[i] = parameter_handler.get_double(axis_labels[i]);
            }
          boundary_parameters.pinned_point_list
            .emplace(index, std::make_pair(parameter_handler.get_double("value"), point));
        }
      else
        {
          AssertThrow(false, FeatureNotImplemented("Vector pinned points"));
        }
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
          std::string subsection_text = "linear solver parameters: ";
          subsection_text.append(variable.name);
          parameter_handler.enter_subsection(subsection_text);

          // Set the tolerance type
          const std::string type_string = parameter_handler.get("tolerance type");
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
          else
            {
              AssertThrow(false, UnreachableCode());
            }

          // Set the tolerance value
          linear_solve_parameters.linear_solve[index].tolerance =
            parameter_handler.get_double("tolerance value");

          // Set the maximum number of iterations
          linear_solve_parameters.linear_solve[index].max_iterations =
            parameter_handler.get_integer("max iterations");

          // Set preconditioner type and related parameters
          linear_solve_parameters.linear_solve[index].preconditioner =
            boost::iequals(parameter_handler.get("preconditioner type"), "GMG")
              ? preconditionerType::GMG
              : preconditionerType::NONE;

          linear_solve_parameters.linear_solve[index].smoothing_range =
            parameter_handler.get_double("smoothing range");

          linear_solve_parameters.linear_solve[index].smoother_degree =
            parameter_handler.get_integer("smoother degree");

          linear_solve_parameters.linear_solve[index].eig_cg_n_iterations =
            parameter_handler.get_integer("eigenvalue cg iterations");

          parameter_handler.leave_subsection();
        }
    }
}

template <int dim>
void
userInputParameters<dim>::assign_nonlinear_solve_parameters(
  dealii::ParameterHandler &parameter_handler)
{
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.field_solve_type == fieldSolveType::NONEXPLICIT_SELF_NONLINEAR ||
          variable.field_solve_type == fieldSolveType::NONEXPLICIT_CO_NONLINEAR)
        {
          std::string subsection_text = "nonlinear solver parameters: ";
          subsection_text.append(variable.name);
          parameter_handler.enter_subsection(subsection_text);

          nonlinear_solve_parameters.nonlinear_solve[index].max_iterations =
            parameter_handler.get_integer("max iterations");
          nonlinear_solve_parameters.nonlinear_solve[index].step_length =
            parameter_handler.get_double("step size");

          // TODO: Implement backtracking line search

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

PRISMS_PF_END_NAMESPACE