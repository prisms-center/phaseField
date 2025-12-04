// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/patterns.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/input_file_reader.h>

#include <prismspf/config.h>

#include <cfloat>
#include <climits>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

InputFileReader::InputFileReader(
  std::string                                       input_file_name,
  const std::map<unsigned int, VariableAttributes> &_var_attributes)
  : parameters_file_name(input_file_name)
  , var_attributes(&_var_attributes)
{
  model_constant_names     = get_model_constant_names();
  const auto num_constants = static_cast<unsigned int>(model_constant_names.size());

  ConditionalOStreams::pout_base()
    << "Number of constants: " << num_constants << "\n"
    << "Number of variables: " << var_attributes->size() << "\n";

  // Read in all of the parameters now
  declare_parameters();
  parameter_handler.parse_input(parameters_file_name);
  number_of_dimensions = static_cast<unsigned int>(parameter_handler.get_integer("dim"));
}

void
InputFileReader::strip_spaces(std::string &line)
{
  while ((!line.empty()) && (line[0] == ' ' || line[0] == '\t'))
    {
      line.erase(0, 1);
    }
  while ((!line.empty()) &&
         (line[line.size() - 1] == ' ' || line[line.size() - 1] == '\t'))
    {
      line.erase(line.size() - 1, std::string::npos);
    }
}

bool
InputFileReader::check_keyword_match(const std::string &line, const std::string &keyword)
{
  // Early return if the line is less than the keyword size
  if (line.size() < keyword.size())
    {
      return false;
    }

  // Check that the line begins with the keyword
  for (unsigned int i = 0; i < keyword.size(); i++)
    {
      if (line[i] != keyword[i])
        {
          return false;
        }
    }

  return true;
}

bool
InputFileReader::parse_line(std::string        line,
                            const std::string &keyword,
                            const std::string &entry_name,
                            std::string       &out_string,
                            const bool         expect_equals_sign)
{
  // Remove spaces from the front and back
  strip_spaces(line);

  // Check whether the line starts with 'keyword'. If not, try next line (if the entry is
  // "", then zero spaces after the keyword is ok)
  if (!check_keyword_match(line, keyword))
    {
      return false;
    }

  if (!entry_name.empty())
    {
      if (line[keyword.size()] != ' ' && line[keyword.size()] != '\t')
        {
          return false;
        }
    }

  // Delete the "keyword" and any more spaces, if present
  line.erase(0, keyword.size());
  strip_spaces(line);

  // now see whether the next word is the word we look for
  if (!line.starts_with(entry_name))
    {
      return false;
    }

  line.erase(0, entry_name.size());
  strip_spaces(line);

  // we'd expect an equals size here if expect_equals_sign is true
  if (expect_equals_sign)
    {
      if ((line.empty()) || (line[0] != '='))
        {
          return false;
        }
    }

  // remove comment
  const std::string::size_type pos = line.find('#');
  if (pos != std::string::npos)
    {
      line.erase(pos);
    }

  // trim the equals sign at the beginning and possibly following spaces
  // as well as spaces at the end
  if (expect_equals_sign)
    {
      line.erase(0, 1);
    }
  strip_spaces(line);

  out_string = line;
  return true;
}

std::set<std::string>
InputFileReader::get_model_constant_names()
{
  const std::string keyword             = "set";
  const std::string entry_name_begining = "Model constant";
  std::ifstream     input_file;
  input_file.open(parameters_file_name);

  std::string line;
  std::string entry;

  std::set<std::string> entry_name_end_list;

  // Loop through each line
  while (std::getline(input_file, line))
    {
      if (parse_line(line, keyword, entry_name_begining, entry, false))
        {
          // Strip whitespace, the equals sign, and everything after the equals
          // sign

          // Strip whitespace at the beginning
          while ((!entry.empty()) && (entry[0] == ' ' || entry[0] == '\t'))
            {
              entry.erase(0, 1);
            }

          // Strip everything up to the equals sign
          while ((!entry.empty()) && (entry[entry.size() - 1] != '='))
            {
              entry.erase(entry.size() - 1, std::string::npos);
            }

          // Strip the equals sign
          entry.erase(entry.size() - 1, std::string::npos);

          // Strip whitespace between the entry name and the equals sign
          while ((!entry.empty()) &&
                 (entry[entry.size() - 1] == ' ' || entry[entry.size() - 1] == '\t'))
            {
              entry.erase(entry.size() - 1, std::string::npos);
            }

          // Add it to the list
          AssertThrow(entry_name_end_list.insert(entry).second,
                      dealii::ExcMessage(
                        "Non-unique constant name in parameters.prm. The constant that "
                        "you attempted to create was \"" +
                        entry + "\"."));
        }
    }
  return entry_name_end_list;
}

void
InputFileReader::declare_parameters()
{
  // Declare all of the entries
  declare_mesh();
  declare_time_discretization();
  declare_solver_parameters();
  declare_output_parameters();
  declare_load_ic_parameters();
  declare_checkpoint_parameters();
  declare_bc_parameters();
  declare_pinning_parameters();
  declare_nucleation_parameters();
  declare_grain_remapping_parameters();
  declare_grain_loading_parameters();
  declare_model_constants();
  declare_miscellaneous_parameters();
}

void
InputFileReader::declare_mesh()
{
  parameter_handler.declare_entry("dim",
                                  "1",
                                  dealii::Patterns::Integer(1, 3),
                                  "The number of dimensions for the simulation.",
                                  true);
  parameter_handler.declare_entry("degree",
                                  "1",
                                  dealii::Patterns::Integer(1,
                                                            Numbers::max_element_degree),
                                  "The polynomial order of the finite element.",
                                  true);
  parameter_handler.declare_entry("global refinement",
                                  "0",
                                  dealii::Patterns::Integer(0, INT_MAX),
                                  "The number of initial refinements of the coarse mesh.",
                                  true);

  parameter_handler.enter_subsection("Rectangular mesh");
  {
    parameter_handler.declare_entry("x size",
                                    "0.0",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The size of the domain in the x direction.");
    parameter_handler.declare_entry("y size",
                                    "0.0",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The size of the domain in the y direction.");
    parameter_handler.declare_entry("z size",
                                    "0.0",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The size of the domain in the z direction.");
    parameter_handler.declare_entry("x lower bound",
                                    "0.0",
                                    dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                    "The lower bound of the domain in the x direction.");
    parameter_handler.declare_entry("y lower bound",
                                    "0.0",
                                    dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                    "The lower bound of the domain in the y direction.");
    parameter_handler.declare_entry("z lower bound",
                                    "0.0",
                                    dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                    "The lower bound of the domain in the z direction.");
    parameter_handler.declare_entry(
      "x subdivisions",
      "1",
      dealii::Patterns::Integer(1, INT_MAX),
      "The number of mesh subdivisions in the x direction.");
    parameter_handler.declare_entry(
      "y subdivisions",
      "1",
      dealii::Patterns::Integer(1, INT_MAX),
      "The number of mesh subdivisions in the y direction.");
    parameter_handler.declare_entry(
      "z subdivisions",
      "1",
      dealii::Patterns::Integer(1, INT_MAX),
      "The number of mesh subdivisions in the z direction.");
  }
  parameter_handler.leave_subsection();

  parameter_handler.enter_subsection("Spherical mesh");
  {
    parameter_handler.declare_entry("radius",
                                    "0",
                                    dealii::Patterns::Double(0.0, DBL_MAX),
                                    "The radius of the domain.");
  }
  parameter_handler.leave_subsection();

  parameter_handler.declare_entry("mesh adaptivity",
                                  "false",
                                  dealii::Patterns::Bool(),
                                  "Whether to enable mesh adaptivity.");
  parameter_handler.declare_entry("max refinement",
                                  "0",
                                  dealii::Patterns::Integer(0, INT_MAX),
                                  "The maximum level of refinement.");
  parameter_handler.declare_entry("min refinement",
                                  "0",
                                  dealii::Patterns::Integer(0, INT_MAX),
                                  "The minimum level of refinement.");
  parameter_handler.declare_entry(
    "remeshing period",
    "2147483647",
    dealii::Patterns::Integer(1, INT_MAX),
    "The number of time steps between mesh refinement operations.");

  for (const auto &[index, variable] : *var_attributes)
    {
      std::string subsection_text = "refinement criterion: ";
      subsection_text.append(variable.get_name());
      parameter_handler.enter_subsection(subsection_text);
      {
        parameter_handler.declare_entry(
          "type",
          "none",
          dealii::Patterns::Selection("none|value|gradient|value_and_gradient"),
          "The type of criterion used to determine if a cell should be "
          "refined. The options are none, value, gradient, value_and_gradient.");
        parameter_handler.declare_entry(
          "value lower bound",
          "0.0",
          dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
          "The lower bound for the window determining where the mesh should be "
          "refined.");
        parameter_handler.declare_entry(
          "value upper bound",
          "0.0",
          dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
          "The upper bound for the window determining where the mesh should be "
          "refined.");
        parameter_handler.declare_entry("gradient magnitude lower bound",
                                        "2147483647",
                                        dealii::Patterns::Double(0.0, DBL_MAX),
                                        "The magnitude of the gradient above "
                                        "which the mesh should be refined.");
      }
      parameter_handler.leave_subsection();
    }
}

void
InputFileReader::declare_time_discretization()
{
  parameter_handler.declare_entry("number steps",
                                  "1",
                                  dealii::Patterns::Integer(1, INT_MAX),
                                  "The total number of step for the simulation.");
  parameter_handler.declare_entry("time step",
                                  "0.0",
                                  dealii::Patterns::Double(0.0, DBL_MAX),
                                  "The time step size for the simulation.");
  parameter_handler.declare_entry(
    "end time",
    "0.0",
    dealii::Patterns::Double(0.0, DBL_MAX),
    "The value of simulated time where the simulation ends.");
}

void
InputFileReader::declare_solver_parameters()
{
  // For linear solves
  for (const auto &[index, variable] : *var_attributes)
    {
      if (variable.get_pde_type() == PDEType::TimeIndependent ||
          variable.get_pde_type() == PDEType::ImplicitTimeDependent)
        {
          std::string subsection_text = "linear solver parameters: ";
          subsection_text.append(variable.get_name());
          parameter_handler.enter_subsection(subsection_text);
          {
            parameter_handler.declare_entry("tolerance type",
                                            "AbsoluteResidual",
                                            dealii::Patterns::Selection(
                                              "AbsoluteResidual|RelativeResidualChange"),
                                            "The tolerance type for the linear solver.");
            parameter_handler.declare_entry(
              "tolerance value",
              "1.0e-10",
              dealii::Patterns::Double(DBL_MIN, DBL_MAX),
              "The value of for the linear solver tolerance.");
            parameter_handler.declare_entry(
              "max iterations",
              "100",
              dealii::Patterns::Integer(1, INT_MAX),
              "The maximum number of linear solver iterations before the loop "
              "is stopped.");
            parameter_handler.declare_entry(
              "preconditioner type",
              "GMG",
              dealii::Patterns::Selection("None|GMG"),
              "The preconditioner type for the linear solver.");
            parameter_handler.declare_entry("smoothing range",
                                            "15.0",
                                            dealii::Patterns::Double(DBL_MIN, DBL_MAX),
                                            "The smoothing range for eigenvalues.");
            parameter_handler.declare_entry("smoother degree",
                                            "5",
                                            dealii::Patterns::Integer(1, INT_MAX),
                                            "The smoother polynomial degree.");
            parameter_handler.declare_entry(
              "eigenvalue cg iterations",
              "10",
              dealii::Patterns::Integer(1, INT_MAX),
              "The maximum number of CG iterations used to find the maximum eigenvalue.");
            parameter_handler.declare_entry("min mg level",
                                            "0",
                                            dealii::Patterns::Integer(0, INT_MAX),
                                            "The minimum multigrid level.");
          }
          parameter_handler.leave_subsection();
        }
    }

  // For nonlinear solves
  for (const auto &[index, variable] : *var_attributes)
    {
      if (variable.get_field_solve_type() == FieldSolveType::NonexplicitSelfnonlinear ||
          variable.get_field_solve_type() == FieldSolveType::NonexplicitCononlinear)
        {
          std::string subsection_text = "nonlinear solver parameters: ";
          subsection_text.append(variable.get_name());
          parameter_handler.enter_subsection(subsection_text);
          {
            parameter_handler.declare_entry("max iterations",
                                            "100",
                                            dealii::Patterns::Integer(1, INT_MAX),
                                            "The maximum number of nonlinear solver "
                                            "iterations before the loop is stopped.");
            parameter_handler.declare_entry(
              "tolerance type",
              "AbsoluteResidual",
              dealii::Patterns::Selection("AbsoluteResidual|RelativeResidualChange"),
              "The tolerance type for the nonlinear solver.");
            parameter_handler.declare_entry(
              "tolerance value",
              "1.0e-10",
              dealii::Patterns::Double(DBL_MIN, DBL_MAX),
              "The value of for the nonlinear solver tolerance.");
            parameter_handler.declare_entry(
              "use backtracking line search",
              "true",
              dealii::Patterns::Bool(),
              "Whether to use a backtracking line-search to find the best "
              "choice of the damping coefficient.");
            parameter_handler.declare_entry(
              "step size modifier",
              "0.5",
              dealii::Patterns::Double(0.0, 1.0),
              "The constant that determines how much the step size decreases "
              "per backtrack. The 'tau' parameter.");
            parameter_handler.declare_entry(
              "residual decrease coefficient",
              "0.5",
              dealii::Patterns::Double(0.0, 1.0),
              "The constant that determines how much the residual must "
              "decrease to be accepted as sufficient. The 'c' parameter.");
            parameter_handler.declare_entry(
              "step size",
              "1.0",
              dealii::Patterns::Double(0.0, 1.0),
              "The constant damping value to be used if the backtrace "
              "line-search approach isn't used.");
          }
          parameter_handler.leave_subsection();
        }
    }
}

void
InputFileReader::declare_output_parameters()
{
  parameter_handler.enter_subsection("output");
  {
    parameter_handler.declare_entry("file name",
                                    "solution",
                                    dealii::Patterns::Anything(),
                                    "The name for the output file, before the "
                                    "time step and processor info are added.");
    parameter_handler.declare_entry("file type",
                                    "vtu",
                                    dealii::Patterns::Selection("vtu|vtk|pvtu"),
                                    "The output file type (either vtu, pvtu, or vtk).");
    parameter_handler.declare_entry(
      "subdivisions",
      "0",
      dealii::Patterns::Integer(0, INT_MAX),
      "The number of subdivisions to apply to the mesh when building output patches.");
    parameter_handler.declare_entry(
      "condition",
      "EQUAL_SPACING",
      dealii::Patterns::Selection("EQUAL_SPACING|LOG_SPACING|N_PER_DECADE|LIST"),
      "The spacing type for outputing the solution fields (either EQUAL_SPACING, "
      "LOG_SPACING, N_PER_DECADE, or LIST).");
    parameter_handler.declare_entry(
      "list",
      "0",
      dealii::Patterns::List(dealii::Patterns::Integer(0, INT_MAX), 0, INT_MAX, ","),
      "The list of time steps to output, used for the LIST type.");
    parameter_handler.declare_entry("number",
                                    "10",
                                    dealii::Patterns::Integer(0, INT_MAX),
                                    "The number of outputs (or number of outputs "
                                    "per decade for the N_PER_DECADE type).");
    parameter_handler.declare_entry(
      "print step period",
      "2147483647",
      dealii::Patterns::Integer(1, INT_MAX),
      "The number of time steps between updates to the screen.");
    parameter_handler.declare_entry(
      "timing information with output",
      "false",
      dealii::Patterns::Bool(),
      "Whether to print the summary table of the wall time and wall time for "
      "indiviual subroutines every time the code outputs.");
  }
  parameter_handler.leave_subsection();
}

void
InputFileReader::declare_load_ic_parameters()
{
  parameter_handler.declare_entry("read initial conditions from file",
                                  "false",
                                  dealii::Patterns::Bool(),
                                  "Whether to read any initial conditions from file.");

  for (unsigned int i = 0; i < Numbers::max_subsections; i++)
    {
      parameter_handler.enter_subsection("initial condition file " + std::to_string(i));
      {
        parameter_handler.declare_entry("file name",
                                        "",
                                        dealii::Patterns::Anything(),
                                        "The file name to load from for each variable.");
        parameter_handler.declare_entry("dataset format",
                                        "vtk_unstructured_grid",
                                        dealii::Patterns::Anything(),
                                        "The type of grid in the file.");
        parameter_handler.declare_entry("file variable names",
                                        "",
                                        dealii::Patterns::List(
                                          dealii::Patterns::Anything()),
                                        "The name of the variable in the file.");
        parameter_handler.declare_entry("simulation variable names",
                                        "",
                                        dealii::Patterns::List(
                                          dealii::Patterns::Anything()),
                                        "The name of the variable in the file.");
        parameter_handler.declare_entry(
          "data points in x direction",
          "-1",
          dealii::Patterns::Integer(-1, INT_MAX),
          "The number of data points of the input file in the x direction.");
        parameter_handler.declare_entry(
          "data points in y direction",
          "-1",
          dealii::Patterns::Integer(-1, INT_MAX),
          "The number of data points of the input file in the y direction.");
        parameter_handler.declare_entry(
          "data points in z direction",
          "-1",
          dealii::Patterns::Integer(-1, INT_MAX),
          "The number of data points of the input file in the z direction.");
      }
      parameter_handler.leave_subsection();
    }
}

void
InputFileReader::declare_checkpoint_parameters()
{
  parameter_handler.enter_subsection("checkpoints");
  {
    parameter_handler.declare_entry(
      "load from checkpoint",
      "false",
      dealii::Patterns::Bool(),
      "Whether to load from a checkpoint created during a previous simulation.");
    parameter_handler.declare_entry(
      "condition",
      "EQUAL_SPACING",
      dealii::Patterns::Selection("EQUAL_SPACING|LOG_SPACING|N_PER_DECADE|LIST"),
      "The spacing type for saving checkpoints (either EQUAL_SPACING, "
      "LOG_SPACING, N_PER_DECADE, or LIST).");
    parameter_handler.declare_entry(
      "list",
      "0",
      dealii::Patterns::List(dealii::Patterns::Integer(0, INT_MAX), 0, INT_MAX, ","),
      "The list of time steps to save checkpoints, used for the LIST type.");
    parameter_handler.declare_entry(
      "number",
      "0",
      dealii::Patterns::Integer(0, INT_MAX),
      "The number of checkpoints (or number of checkpoints per decade for the "
      "N_PER_DECADE type).");
  }
  parameter_handler.leave_subsection();
}

void
InputFileReader::declare_bc_parameters()
{
  for (const auto &[index, variable] : *var_attributes)
    {
      if (variable.is_postprocess())
        {
          continue;
        }
      if (variable.field_info.tensor_rank == FieldInfo::TensorRank::Scalar)
        {
          std::string bc_text = "boundary condition for ";
          bc_text.append(variable.get_name());
          parameter_handler.declare_entry(
            bc_text,
            "",
            dealii::Patterns::Anything(),
            "The boundary conditions for one of the governing equations).");
        }
      else
        {
          std::string bc_text = "boundary condition for ";
          bc_text.append(variable.get_name());
          bc_text.append(", x component");
          parameter_handler.declare_entry(
            bc_text,
            "",
            dealii::Patterns::Anything(),
            "The boundary conditions for one of the governing equations).");

          bc_text = "boundary condition for ";
          bc_text.append(variable.get_name());
          bc_text.append(", y component");
          parameter_handler.declare_entry(
            bc_text,
            "",
            dealii::Patterns::Anything(),
            "The boundary conditions for one of the governing equations).");

          bc_text = "boundary condition for ";
          bc_text.append(variable.get_name());
          bc_text.append(", z component");
          parameter_handler.declare_entry(
            bc_text,
            "",
            dealii::Patterns::Anything(),
            "The boundary conditions for one of the governing equations).");
        }
    }
}

void
InputFileReader::declare_pinning_parameters()
{
  for (const auto &[index, variable] : *var_attributes)
    {
      if (variable.is_postprocess())
        {
          continue;
        }
      std::string pinning_text = "pinning point for ";
      pinning_text.append(variable.get_name());
      parameter_handler.enter_subsection(pinning_text);
      {
        if (variable.field_info.tensor_rank == FieldInfo::TensorRank::Scalar)
          {
            parameter_handler.declare_entry("value",
                                            "2147483647",
                                            dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                            "Value of pinned point.");
          }
        else
          {
            parameter_handler.declare_entry("x value",
                                            "2147483647",
                                            dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                            "Value of pinned point for the x-component.");
            parameter_handler.declare_entry("y value",
                                            "2147483647",
                                            dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                            "Value of pinned point for the y-component.");
            parameter_handler.declare_entry("z value",
                                            "2147483647",
                                            dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                            "Value of pinned point for the z-component.");
          }
        parameter_handler.declare_entry("x",
                                        "0.0",
                                        dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                        "X-coordinate of the point");
        parameter_handler.declare_entry("y",
                                        "0.0",
                                        dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                        "Y-coordinate of the point");
        parameter_handler.declare_entry("z",
                                        "0.0",
                                        dealii::Patterns::Double(-DBL_MAX, DBL_MAX),
                                        "Z-coordinate of the point");
      }
      parameter_handler.leave_subsection();
    }
}

void
InputFileReader::declare_nucleation_parameters()
{
  parameter_handler.enter_subsection("nucleation");
  {
    parameter_handler.declare_entry("nucleus exclusion distance",
                                    "0.0",
                                    dealii::Patterns::Double(),
                                    "The minimum distance between nuclei.");
    parameter_handler.declare_entry("same field nucleus exclusion distance",
                                    "0.0",
                                    dealii::Patterns::Double(),
                                    "The minimum distance between nuclei.");
    parameter_handler.declare_entry(
      "nucleation period",
      "1",
      dealii::Patterns::Integer(1),
      "The number of increments between nucleation attempts.");
    parameter_handler.declare_entry(
      "refinement radius",
      "0.0",
      dealii::Patterns::Double(0.0),
      "The radius around a nucleus in which AMR is applied.");
    parameter_handler.declare_entry(
      "seeding time",
      "0.0",
      dealii::Patterns::Double(0.0),
      "The time duration over which nuclei are considered \"active\" and refinement and "
      "exclusion zones are applied. Same as \"seeding increments\" but in time.");
    parameter_handler.declare_entry(
      "seeding increments",
      "1",
      dealii::Patterns::Integer(1, INT_MAX),
      "The number of increments over which nuclei are considered \"active\" and "
      "refinement and exclusion zones are applied. Same as \"seeding time\" but in "
      "increments.");
    { // Declare aliases for the parameters
      //============================================================================================
      parameter_handler.declare_alias("nucleus exclusion distance",
                                      "nucleus_exclusion_distance");
      parameter_handler.declare_alias("nucleus exclusion distance",
                                      "nucleus exclusion radius");
      parameter_handler.declare_alias("nucleus exclusion distance",
                                      "nucleus_exclusion_radius");
      parameter_handler.declare_alias("nucleus exclusion distance", "exclusion distance");
      parameter_handler.declare_alias("nucleus exclusion distance", "exclusion_distance");
      parameter_handler.declare_alias("nucleus exclusion distance", "exclusion radius");
      parameter_handler.declare_alias("nucleus exclusion distance", "exclusion_radius");
      //
      parameter_handler.declare_alias("same field nucleus exclusion distance",
                                      "same_field_nucleus_exclusion_distance");
      parameter_handler.declare_alias("same field nucleus exclusion distance",
                                      "same field nucleus exclusion radius");
      parameter_handler.declare_alias("same field nucleus exclusion distance",
                                      "same_field_nucleus_exclusion_radius");
      parameter_handler.declare_alias("same field nucleus exclusion distance",
                                      "same field exclusion distance");
      parameter_handler.declare_alias("same field nucleus exclusion distance",
                                      "same_field_exclusion_distance");
      parameter_handler.declare_alias("same field nucleus exclusion distance",
                                      "same field exclusion radius");
      parameter_handler.declare_alias("same field nucleus exclusion distance",
                                      "same_field_exclusion_radius");
    }
  }
  parameter_handler.leave_subsection();
}

void
InputFileReader::declare_grain_remapping_parameters()
{
  // parameter_handler.declare_entry(
  //   "Activate grain reassignment",
  //   "false",
  //   dealii::Patterns::Bool(),
  //   "Whether to enable the grain reassignment capabilities of PRISMS-PF where "
  //   "multiple grains are packed into a single order parameter.");

  // parameter_handler.declare_entry(
  //   "Time steps between grain reassignments",
  //   "100",
  //   dealii::Patterns::Integer(),
  //   "The number of time steps between times when the grain reassignment "
  //   "algorithm is triggered.");

  // parameter_handler.declare_entry(
  //   "Order parameter cutoff for grain identification",
  //   "1.0e-4",
  //   dealii::Patterns::Double(),
  //   "The threshold value of the order parameter where the element is "
  //   "considered to be in the grain or out of the grain.");

  // parameter_handler.declare_entry(
  //   "Buffer between grains before reassignment",
  //   "-1.0",
  //   dealii::Patterns::Double(),
  //   "The buffer value added to the radius of all grains used to calculation "
  //   "whether grains should be reassigned.");

  // parameter_handler.declare_entry(
  //   "Order parameter fields for grain reassignment",
  //   "",
  //   dealii::Patterns::List(dealii::Patterns::Anything()),
  //   "The list of field indices for the shared order parameters for grain "
  //   "reassignment.");
}

void
InputFileReader::declare_grain_loading_parameters()
{
  // parameter_handler.declare_entry("Load grain structure",
  //                                 "false",
  //                                 dealii::Patterns::Bool(),
  //                                 "Whether to load a grain structure in from file.");

  // parameter_handler.declare_entry(
  //   "vtk file type",
  //   "UNSTRUCTURED",
  //   dealii::Patterns::Anything(),
  //   "Whether to load an unstructured file for grain structure."); // reads the type
  //   of
  //                                                                 // file from the
  //                                                                 input
  //                                                                 // parameters.prm
  //                                                                 file,
  //                                                                 // deafault setting
  //                                                                 is
  //                                                                 // unstructured
  //                                                                 mesh

  // parameter_handler.declare_entry(
  //   "Grain structure filename",
  //   "",
  //   dealii::Patterns::Anything(),
  //   "The filename (not including the '.vtk' extension) for the file holding "
  //   "the grain structure to be loaded.");

  // parameter_handler.declare_entry(
  //   "Grain structure variable name",
  //   "",
  //   dealii::Patterns::Anything(),
  //   "The variable name in the file holding the grain structure to be loaded "
  //   "that contains the grain ids.");

  // parameter_handler.declare_entry(
  //   "Number of smoothing cycles after grain structure loading",
  //   "10",
  //   dealii::Patterns::Integer(),
  //   "The number of times a diffusion smoother is run on the order parameters "
  //   "after the grains are loaded from file. The smoothing is necessary for the "
  //   "adaptive mesher to work properly.");

  // parameter_handler.declare_entry(
  //   "Minimum radius for loaded grains",
  //   "0.0",
  //   dealii::Patterns::Double(),
  //   "The minimum radius for a body to be considered a grain instead of an "
  //   "artifact from the loading process.");
}

void
InputFileReader::declare_miscellaneous_parameters()
{
  parameter_handler.enter_subsection("miscellaneous");
  {
    parameter_handler.declare_entry(
      "random seed",
      "2025",
      dealii::Patterns::Integer(0, INT_MAX),
      "The random seed for the simulation. "
      "This is used to initialize the random number generator.");
  }
  parameter_handler.leave_subsection();
}

void
InputFileReader::declare_model_constants()
{
  for (const std::string &constant_name : model_constant_names)
    {
      std::string constants_text = "Model constant ";
      constants_text.append(constant_name);
      parameter_handler.declare_entry(constants_text,
                                      "0",
                                      dealii::Patterns::Anything(),
                                      "The value of a user-defined constant.");
    }
}

PRISMS_PF_END_NAMESPACE
