#include <deal.II/base/mpi.h>
#include <deal.II/base/utilities.h>

#include <core/exceptions.h>
#include <core/inputFileReader.h>
#include <core/refinement/RefinementCriterion.h>
#include <fstream>

inputFileReader::inputFileReader(const std::string    &input_file_name,
                                 const AttributesList &_var_attributes,
                                 const AttributesList &_pp_attributes)
  : parameters_file_name(input_file_name)
  , var_attributes(_var_attributes)
  , pp_attributes(_pp_attributes)
{
  model_constant_names = get_model_constant_names();
  uint num_constants   = model_constant_names.size();

  if (dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      std::cout << "Number of constants: " << num_constants << "\n";
      std::cout << "Number of post-processing variables: " << pp_attributes.size()
                << "\n";
    }

  // Read in all of the parameters now
  declare_parameters();
  parameter_handler.parse_input(parameters_file_name);
  number_of_dimensions = parameter_handler.get_integer("Number of dimensions");
}

void
inputFileReader::strip_spaces(std::string &line)
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
inputFileReader::check_keyword_match(std::string &line, const std::string &keyword)
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
inputFileReader::parse_line(std::string        line,
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
  if (line.find(entry_name) != 0)
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

// Method to parse an input file to get a list of variables from related
// subsections
std::vector<std::string>
inputFileReader::get_subsection_entry_list(const std::string &subsec_name,
                                           const std::string &entry_name,
                                           const std::string &default_entry)
{
  std::ifstream input_file;
  input_file.open(parameters_file_name);

  std::string               line;
  std::string               entry;
  bool                      in_subsection       = false;
  bool                      found_entry         = false;
  bool                      desired_entry_found = false;
  unsigned int              subsection_index    = 0;
  std::vector<std::string>  entry_list;
  std::vector<unsigned int> index_list;

  // Loop through each line
  while (std::getline(input_file, line))
    {
      // If the line is the start of a subsection, turn 'in_subsection' to true
      // and store the subsection index
      if (!in_subsection)
        {
          found_entry = parse_line(line, "subsection", subsec_name, entry, false);
          if (found_entry)
            {
              in_subsection       = true;
              subsection_index    = dealii::Utilities::string_to_int(entry);
              desired_entry_found = false;
            }
        }
      // If in a subsection, look for the line setting the entry or for the end
      // of the subsection
      else
        {
          found_entry = parse_line(line, "set", entry_name, entry, true);
          if (found_entry)
            {
              entry_list.push_back(entry);
              index_list.push_back(subsection_index);
              desired_entry_found = true;
            }
          found_entry = parse_line(line, "end", "", entry, false);
          if (found_entry)
            {
              if (!desired_entry_found)
                {
                  entry_list.push_back(default_entry);
                  index_list.push_back(subsection_index);
                }
              in_subsection       = false;
              desired_entry_found = false;
            }
        }
    }

  // Now sort the entry list vector so that it is in index order
  std::vector<std::string> sorted_entry_list;
  for (unsigned int i = 0; i < entry_list.size(); i++)
    {
      for (unsigned int j = 0; j < entry_list.size(); j++)
        {
          if (i == j)
            {
              sorted_entry_list.push_back(entry_list[index_list[j]]);
              break;
            }
        }
    }
  return sorted_entry_list;
}

// Method to parse an input file to get a list of variables from related
// subsections
std::set<std::string>
inputFileReader::get_model_constant_names()
{
  const std::string keyword             = "set";
  const std::string entry_name_begining = "Model constant";
  std::ifstream     input_file;
  input_file.open(parameters_file_name);

  std::string line;
  std::string entry;
  bool        found_entry = false;

  std::set<std::string> entry_name_end_list;

  // Loop through each line
  while (std::getline(input_file, line))
    {
      found_entry = parse_line(line, keyword, entry_name_begining, entry, false);
      if (found_entry)
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
          AssertThrow(
            entry_name_end_list.insert(entry).second,
            dealii::ExcMessage(
              "PRISMS-PF Error: Non-unique constant name in parameters.prm. The "
              "constant that you attempted to create was \"" +
              entry + "\"."));
        }
    }
  return entry_name_end_list;
}

void
inputFileReader::declare_parameters()
{
  // Declare all of the entries
  declare_mesh();
  declare_time_discretization();
  declare_solver_parameters();
  declare_output_parameters();
  declare_load_IC_parameters();
  declare_checkpoint_parameters();
  declare_BC_parameters();
  declare_pinning_parameters();
  declare_nucleation_parameters();
  declare_grain_remapping_parameters();
  declare_grain_loading_parameters();
  declare_model_constants();
}

void
inputFileReader::declare_mesh()
{
  parameter_handler.declare_entry("Number of dimensions",
                                  "-1",
                                  dealii::Patterns::Integer(),
                                  "The number of dimensions for the simulation.");

  parameter_handler.declare_entry("Domain size X",
                                  "-1",
                                  dealii::Patterns::Double(),
                                  "The size of the domain in the x direction.");
  parameter_handler.declare_entry("Domain size Y",
                                  "-1",
                                  dealii::Patterns::Double(),
                                  "The size of the domain in the y direction.");
  parameter_handler.declare_entry("Domain size Z",
                                  "-1",
                                  dealii::Patterns::Double(),
                                  "The size of the domain in the z direction.");
  parameter_handler.declare_entry("Element degree",
                                  "1",
                                  dealii::Patterns::Integer(),
                                  "The polynomial order of the finte element.");
  parameter_handler.declare_entry("Subdivisions X",
                                  "1",
                                  dealii::Patterns::Integer(),
                                  "The number of mesh subdivisions in the x direction.");
  parameter_handler.declare_entry("Subdivisions Y",
                                  "1",
                                  dealii::Patterns::Integer(),
                                  "The number of mesh subdivisions in the y direction.");
  parameter_handler.declare_entry("Subdivisions Z",
                                  "1",
                                  dealii::Patterns::Integer(),
                                  "The number of mesh subdivisions in the z direction.");
  parameter_handler.declare_entry(
    "Refine factor",
    "-1",
    dealii::Patterns::Integer(),
    "The number of initial refinements of the coarse mesh.");

  parameter_handler.declare_entry("Mesh adaptivity",
                                  "false",
                                  dealii::Patterns::Bool(),
                                  "Whether to enable mesh adaptivity.");
  parameter_handler.declare_entry("Max refinement level",
                                  "-1",
                                  dealii::Patterns::Integer(),
                                  "The maximum level of refinement.");
  parameter_handler.declare_entry("Min refinement level",
                                  "-1",
                                  dealii::Patterns::Integer(),
                                  "The minimum level of refinement.");
  parameter_handler.declare_entry(
    "Steps between remeshing operations",
    "1",
    dealii::Patterns::Integer(),
    "The number of time steps between mesh refinement operations.");

  for (const auto &[index, variable] : var_attributes)
    {
      std::string subsection_text = "Refinement criterion: ";
      subsection_text.append(variable.name);
      parameter_handler.enter_subsection(subsection_text);
      {
        parameter_handler.declare_entry(
          "Criterion type",
          "",
          dealii::Patterns::Anything(),
          "The type of criterion used to determine if a cell should be "
          "refined. The options are VALUE, GRADIENT, VALUE_AND_GRADIENT.");

        parameter_handler.declare_entry(
          "Value lower bound",
          "0.0",
          dealii::Patterns::Double(),
          "The lower bound for the window determining where the mesh should be "
          "refined.");

        parameter_handler.declare_entry(
          "Value upper bound",
          "1.0",
          dealii::Patterns::Double(),
          "The upper bound for the window determining where the mesh should be "
          "refined.");

        parameter_handler.declare_entry("Gradient magnitude lower bound",
                                        "1.0",
                                        dealii::Patterns::Double(),
                                        "The magnitude of the gradient above "
                                        "which the mesh should be refined.");
      }
      parameter_handler.leave_subsection();
    }
}

void
inputFileReader::declare_time_discretization()
{
  parameter_handler.declare_entry("Number of time steps",
                                  "-1",
                                  dealii::Patterns::Integer(),
                                  "The time step size for the simulation.");
  parameter_handler.declare_entry("Time step",
                                  "-0.1",
                                  dealii::Patterns::Double(),
                                  "The time step size for the simulation.");
  parameter_handler.declare_entry(
    "Simulation end time",
    "-0.1",
    dealii::Patterns::Double(),
    "The value of simulated time where the simulation ends.");
}

void
inputFileReader::declare_solver_parameters()
{
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.eq_type == TIME_INDEPENDENT ||
          variable.eq_type == IMPLICIT_TIME_DEPENDENT)
        {
          std::string subsection_text = "Linear solver parameters: ";
          subsection_text.append(variable.name);
          parameter_handler.enter_subsection(subsection_text);
          {
            parameter_handler.declare_entry("Tolerance type",
                                            "RELATIVE_RESIDUAL_CHANGE",
                                            dealii::Patterns::Anything(),
                                            "The tolerance type for the linear solver.");
            parameter_handler.declare_entry(
              "Tolerance value",
              "1.0e-10",
              dealii::Patterns::Double(),
              "The value of for the linear solver tolerance.");
            parameter_handler.declare_entry(
              "Maximum linear solver iterations",
              "1000",
              dealii::Patterns::Integer(),
              "The maximum number of linear solver iterations before the loop "
              "is stopped.");
          }
          parameter_handler.leave_subsection();
        }
    }

  parameter_handler.declare_entry("Maximum nonlinear solver iterations",
                                  "100",
                                  dealii::Patterns::Integer(),
                                  "The maximum number of nonlinear solver "
                                  "iterations before the loop is stopped.");

  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.is_nonlinear)
        {
          std::string subsection_text = "Nonlinear solver parameters: ";
          subsection_text.append(variable.name);
          parameter_handler.enter_subsection(subsection_text);
          {
            parameter_handler.declare_entry(
              "Tolerance type",
              "ABSOLUTE_SOLUTION_CHANGE",
              dealii::Patterns::Anything(),
              "The tolerance type for the nonlinear solver.");
            parameter_handler.declare_entry(
              "Tolerance value",
              "1.0e-10",
              dealii::Patterns::Double(),
              "The value of for the nonlinear solver tolerance.");
            parameter_handler.declare_entry(
              "Use backtracking line search damping",
              "true",
              dealii::Patterns::Bool(),
              "Whether to use a backtracking line-search to find the best "
              "choice of the damping coefficient.");
            parameter_handler.declare_entry(
              "Backtracking step size modifier",
              "0.5",
              dealii::Patterns::Double(),
              "The constant that determines how much the step size decreases "
              "per backtrack. The 'tau' parameter.");
            parameter_handler.declare_entry(
              "Backtracking residual decrease coefficient",
              "1.0",
              dealii::Patterns::Double(),
              "The constant that determines how much the residual must "
              "decrease to be accepted as sufficient. The 'c' parameter.");
            parameter_handler.declare_entry(
              "Constant damping value",
              "1.0",
              dealii::Patterns::Double(),
              "The constant damping value to be used if the backtrace "
              "line-search approach isn't used.");
            parameter_handler.declare_entry(
              "Use Laplace's equation to determine the initial guess",
              "false",
              dealii::Patterns::Bool(),
              "Whether to use the solution of Laplace's equation instead of "
              "the IC in ICs_and_BCs.h as the initial guess for nonlinear, "
              "time independent equations. This guarantees smoothness and "
              "compliance with BCs.");
          }
          parameter_handler.leave_subsection();
        }
    }
}

void
inputFileReader::declare_output_parameters()
{
  parameter_handler.declare_entry("Output file name (base)",
                                  "solution",
                                  dealii::Patterns::Anything(),
                                  "The name for the output file, before the "
                                  "time step and processor info are added.");
  parameter_handler.declare_entry("Output file type",
                                  "vtu",
                                  dealii::Patterns::Anything(),
                                  "The output file type (either vtu or vtk).");
  parameter_handler.declare_entry(
    "Output separate files per process",
    "false",
    dealii::Patterns::Bool(),
    "Whether to output separate vtu files for each process in a parallel "
    "calculation (automatically set to true for vtk files).");
  parameter_handler.declare_entry("Output condition",
                                  "EQUAL_SPACING",
                                  dealii::Patterns::Anything(),
                                  "The spacing type for outputing the solution fields.");
  parameter_handler.declare_entry(
    "List of time steps to output",
    "0",
    dealii::Patterns::Anything(),
    "The list of time steps to output, used for the LIST type.");
  parameter_handler.declare_entry("Number of outputs",
                                  "10",
                                  dealii::Patterns::Integer(),
                                  "The number of outputs (or number of outputs "
                                  "per decade for the N_PER_DECADE type).");
  parameter_handler.declare_entry(
    "Skip print steps",
    "1",
    dealii::Patterns::Integer(),
    "The number of time steps between updates to the screen.");
  parameter_handler.declare_entry(
    "Print timing information with output",
    "false",
    dealii::Patterns::Bool(),
    "Whether to print the summary table of the wall time and wall time for "
    "indiviual subroutines every time the code outputs.");
}

void
inputFileReader::declare_load_IC_parameters()
{
  parameter_handler.declare_entry(
    "Load initial conditions",
    "void",
    dealii::Patterns::Anything(),
    "Whether to load the initial conditions for each variable from file.");
  parameter_handler.declare_entry(
    "Load parallel file",
    "void",
    dealii::Patterns::Anything(),
    "Whether all processors should read from a single file (versus each "
    "reading from separate files).");
  parameter_handler.declare_entry("File names",
                                  "void",
                                  dealii::Patterns::Anything(),
                                  "The file name to load from for each variable.");
  parameter_handler.declare_entry(
    "Variable names in the files",
    "void",
    dealii::Patterns::Anything(),
    "What each variable is named in the file being loaded.");
}

void
inputFileReader::declare_checkpoint_parameters()
{
  parameter_handler.declare_entry(
    "Load from a checkpoint",
    "false",
    dealii::Patterns::Bool(),
    "Whether to load from a checkpoint created during a previous simulation.");
  parameter_handler.declare_entry("Checkpoint condition",
                                  "EQUAL_SPACING",
                                  dealii::Patterns::Anything(),
                                  "The spacing type for saving checkpoints.");
  parameter_handler.declare_entry(
    "List of time steps to save checkpoints",
    "0",
    dealii::Patterns::Anything(),
    "The list of time steps to save checkpoints, used for the LIST type.");
  parameter_handler.declare_entry(
    "Number of checkpoints",
    "0",
    dealii::Patterns::Integer(),
    "The number of checkpoints (or number of checkpoints per decade for the "
    "N_PER_DECADE type).");
}

void
inputFileReader::declare_BC_parameters()
{
  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.var_type == SCALAR)
        {
          std::string bc_text = "Boundary condition for variable ";
          bc_text.append(variable.name);
          parameter_handler.declare_entry(
            bc_text,
            "",
            dealii::Patterns::Anything(),
            "The boundary conditions for one of the governing equations).");
        }
      else
        {
          std::string bc_text = "Boundary condition for variable ";
          bc_text.append(variable.name);
          bc_text.append(", x component");
          parameter_handler.declare_entry(
            bc_text,
            "",
            dealii::Patterns::Anything(),
            "The boundary conditions for one of the governing equations).");

          bc_text = "Boundary condition for variable ";
          bc_text.append(variable.name);
          bc_text.append(", y component");
          parameter_handler.declare_entry(
            bc_text,
            "",
            dealii::Patterns::Anything(),
            "The boundary conditions for one of the governing equations).");

          bc_text = "Boundary condition for variable ";
          bc_text.append(variable.name);
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
inputFileReader::declare_pinning_parameters()
{
  for (const auto &[index, variable] : var_attributes)
    {
      std::string pinning_text = "Pinning point: ";
      pinning_text.append(variable.name);
      parameter_handler.enter_subsection(pinning_text);
      {
        parameter_handler.declare_entry("x",
                                        "-1.0",
                                        dealii::Patterns::Double(),
                                        "X-coordinate of the point");
        parameter_handler.declare_entry("y",
                                        "0.0",
                                        dealii::Patterns::Double(),
                                        "Y-coordinate of the point");
        parameter_handler.declare_entry("z",
                                        "0.0",
                                        dealii::Patterns::Double(),
                                        "Z-coordinate of the point");
      }
      parameter_handler.leave_subsection();
    }
}

void
inputFileReader::declare_nucleation_parameters()
{
  parameter_handler.declare_entry(
    "Enable evolution before nucleation",
    "false",
    dealii::Patterns::Bool(),
    "Whether variable fields evolve before the first nucleus appears.");
  // Implement later
  // parameter_handler.declare_entry("Allow multiple nuclei per order
  // parameter","true",dealii::Patterns::Bool(),"Whether multiple nucleation
  // events can occur within an order parameter.");
  parameter_handler.declare_entry("Minimum allowed distance between nuclei",
                                  "-1",
                                  dealii::Patterns::Double(),
                                  "The minimum allowed distance between nuclei "
                                  "placed during the same time step.");
  parameter_handler.declare_entry(
    "Order parameter cutoff value",
    "0.01",
    dealii::Patterns::Double(),
    "Order parameter cutoff value for nucleation (when the sum of all order "
    "parameters is above this value, no nucleation is attempted).");
  parameter_handler.declare_entry(
    "Time steps between nucleation attempts",
    "100",
    dealii::Patterns::Integer(),
    "The number of time steps between nucleation attempts.");
  parameter_handler.declare_entry("Nucleation start time",
                                  "0.0",
                                  dealii::Patterns::Double(),
                                  "The time at which nucleation starts.");
  parameter_handler.declare_entry("Nucleation end time",
                                  "1.0e10",
                                  dealii::Patterns::Double(),
                                  "The time after which no nucleation occurs.");

  for (const auto &[index, variable] : var_attributes)
    {
      if (variable.nucleating_variable)
        {
          std::string nucleation_text = "Nucleation parameters: ";
          nucleation_text.append(variable.name);
          parameter_handler.enter_subsection(nucleation_text);
          {
            parameter_handler.declare_entry(
              "Nucleus semiaxes (x, y, z)",
              "0,0,0",
              dealii::Patterns::List(dealii::Patterns::Double()),
              "The semiaxes for nuclei placed with the explicit nucleation "
              "algorithm.");
            parameter_handler.declare_entry(
              "Nucleus rotation in degrees (x, y, z)",
              "0,0,0",
              dealii::Patterns::List(dealii::Patterns::Double()),
              "The rotation of the nuclei placed with the explicit nucleation "
              "algorithm. The rotations are given with respect to the normal "
              "direction using intrinsic Tait-Bryan angles.");
            parameter_handler.declare_entry(
              "Freeze zone semiaxes (x, y, z)",
              "0,0,0",
              dealii::Patterns::List(dealii::Patterns::Double()),
              "The semiaxes for region where the order parameter is frozen for "
              "a period of time after placement.");
            parameter_handler.declare_entry(
              "Freeze time following nucleation",
              "0.0",
              dealii::Patterns::Double(),
              "Duration that the order parameter is frozen after placement.");
            parameter_handler.declare_entry(
              "Nucleation-free border thickness",
              "0.0",
              dealii::Patterns::Double(),
              "The thickness of the nucleation-free region near the domain "
              "boundaries (ignored for periodic BCs).");
          }
          parameter_handler.leave_subsection();
        }
    }
}

void
inputFileReader::declare_grain_remapping_parameters()
{
  parameter_handler.declare_entry(
    "Activate grain reassignment",
    "false",
    dealii::Patterns::Bool(),
    "Whether to enable the grain reassignment capabilities of PRISMS-PF where "
    "multiple grains are packed into a single order parameter.");

  parameter_handler.declare_entry(
    "Time steps between grain reassignments",
    "100",
    dealii::Patterns::Integer(),
    "The number of time steps between times when the grain reassignment "
    "algorithm is triggered.");

  parameter_handler.declare_entry(
    "Order parameter cutoff for grain identification",
    "1.0e-4",
    dealii::Patterns::Double(),
    "The threshold value of the order parameter where the element is "
    "considered to be in the grain or out of the grain.");

  parameter_handler.declare_entry(
    "Buffer between grains before reassignment",
    "-1.0",
    dealii::Patterns::Double(),
    "The buffer value added to the radius of all grains used to calculation "
    "whether grains should be reassigned.");

  parameter_handler.declare_entry(
    "Order parameter fields for grain reassignment",
    "",
    dealii::Patterns::List(dealii::Patterns::Anything()),
    "The list of field indices for the shared order parameters for grain "
    "reassignment.");
}

void
inputFileReader::declare_grain_loading_parameters()
{
  parameter_handler.declare_entry("Load grain structure",
                                  "false",
                                  dealii::Patterns::Bool(),
                                  "Whether to load a grain structure in from file.");

  parameter_handler.declare_entry(
    "vtk file type",
    "UNSTRUCTURED",
    dealii::Patterns::Anything(),
    "Whether to load an unstructured file for grain structure."); // reads the type of
                                                                  // file from the input
                                                                  // parameters.prm file,
                                                                  // deafault setting is
                                                                  // unstructured mesh

  parameter_handler.declare_entry(
    "Grain structure filename",
    "",
    dealii::Patterns::Anything(),
    "The filename (not including the '.vtk' extension) for the file holding "
    "the grain structure to be loaded.");

  parameter_handler.declare_entry(
    "Grain structure variable name",
    "",
    dealii::Patterns::Anything(),
    "The variable name in the file holding the grain structure to be loaded "
    "that contains the grain ids.");

  parameter_handler.declare_entry(
    "Number of smoothing cycles after grain structure loading",
    "10",
    dealii::Patterns::Integer(),
    "The number of times a diffusion smoother is run on the order parameters "
    "after the grains are loaded from file. The smoothing is necessary for the "
    "adaptive mesher to work properly.");

  parameter_handler.declare_entry(
    "Minimum radius for loaded grains",
    "0.0",
    dealii::Patterns::Double(),
    "The minimum radius for a body to be considered a grain instead of an "
    "artifact from the loading process.");
}

void
inputFileReader::declare_model_constants()
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
