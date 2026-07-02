#include <deal.II/base/exceptions.h>

#include <prismspf/user_inputs/io_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
FieldOutputParameters::predeclare(dealii::ParameterHandler &parameter_handler) const
{
  AssertThrow(false, dealii::ExcNotImplemented());
}

void
FieldOutputParameters::preassign(dealii::ParameterHandler &parameter_handler)
{
  AssertThrow(false, dealii::ExcNotImplemented());
}

void
FieldOutputParameters::declare(dealii::ParameterHandler &parameter_handler,
                               unsigned int              max_criteria) const
{
  parameter_handler.enter_subsection("output");
  {
    parameter_handler.declare_entry(
      "file type",
      "vtu",
      dealii::Patterns::Selection("vtu|vtk|pvtu|xdmf"),
      "The output file type (either vtu, pvtu, vtk, or xdmf).");

    parameter_handler.declare_entry(
      "subdivisions",
      "0",
      dealii::Patterns::Integer(0, INT_MAX),
      "The number of subdivisions to apply to the mesh when building output patches. "
      "If 0, the degree is used.");

    parameter_handler.declare_entry("compression level",
                                    "default",
                                    dealii::Patterns::Selection(
                                      "default|speed|size|none"),
                                    "The compression level for output.");

    parameter_handler.declare_entry("directory",
                                    "outputs",
                                    dealii::Patterns::Anything(),
                                    "The name of the output directory.");
    parameter_handler.declare_alias("directory", "folder name");

    parameter_handler.declare_entry("file name",
                                    "solution",
                                    dealii::Patterns::Anything(),
                                    "The prefix of the output files, before the "
                                    "time step and processor info are added.");

    parameter_handler.declare_entry(
      "condition",
      "EQUAL_SPACING",
      dealii::Patterns::Selection("EQUAL_SPACING|LOG_SPACING|N_PER_DECADE|LIST"),
      "The spacing type for outputting the solution fields.");
    parameter_handler.declare_entry(
      "list",
      "0",
      dealii::Patterns::List(dealii::Patterns::Integer(0, INT_MAX), 0, INT_MAX, ","),
      "The list of time steps to output. Used for the LIST type only and must be comma "
      "delimited.");
    parameter_handler.declare_entry("number",
                                    "10",
                                    dealii::Patterns::Integer(0, INT_MAX),
                                    "The number of outputs (or number of outputs "
                                    "per decade for the N_PER_DECADE type).");

    parameter_handler.declare_entry(
      "variables",
      "",
      dealii::Patterns::List(dealii::Patterns::Anything(), 0, INT_MAX, ","),
      "The list of the fields to output. Must be comma delimited. Additionally, for "
      "the output of left-hand side and old fields, they must follow the same "
      "delimiters that are used in dependency sets. In other words, something like "
      "`set variables = n1, old_1(n1), lhs(n1)`.");
  }
  parameter_handler.leave_subsection();
}

void
FieldOutputParameters::assign(dealii::ParameterHandler &parameter_handler,
                              unsigned int              max_criteria)
{
  AssertThrow(false, dealii::ExcNotImplemented());
}

void
FieldOutputParameters::assign(dealii::ParameterHandler &parameter_handler,
                              unsigned int              n_increments,
                              unsigned int              max_criteria)
{
  const static std::unordered_map<std::string, FieldOutputParameters::OutputType>
    output_type_table = {
      {"vtu",  FieldOutputParameters::OutputType::VTU },
      {"vtk",  FieldOutputParameters::OutputType::VTK },
      {"pvtu", FieldOutputParameters::OutputType::PVTU},
      {"xdmf", FieldOutputParameters::OutputType::XDMF}
  };
  const static std::unordered_map<std::string, dealii::DataOutBase::CompressionLevel>
    compression_level_table = {
      {"default",    dealii::DataOutBase::CompressionLevel::default_compression},
      {"best speed", dealii::DataOutBase::CompressionLevel::best_speed         },
      {"best size",  dealii::DataOutBase::CompressionLevel::best_compression   },
      {"none",       dealii::DataOutBase::CompressionLevel::no_compression     }
  };

  parameter_handler.enter_subsection("output");
  {
    folder             = parameter_handler.get("directory");
    file_name          = parameter_handler.get("file name");
    patch_subdivisions = (unsigned int) (parameter_handler.get_integer("subdivisions"));

    file_type = output_type_table.at(parameter_handler.get("file type"));

    compression_level =
      compression_level_table.at(parameter_handler.get("compression level"));

    add_list_outputs(dealii::Utilities::split_string_list(
                       parameter_handler.get("variables")),
                     output_fields);

    std::string  condition = parameter_handler.get("condition");
    unsigned int n_outputs = (unsigned int) (parameter_handler.get_integer("number"));

    if (condition == "EQUAL_SPACING")
      {
        add_equal_spacing_outputs(n_outputs, n_increments, output_list);
      }
    else if (condition == "LOG_SPACING")
      {
        add_log_spacing_outputs(n_outputs, n_increments, output_list);
      }
    else if (condition == "N_PER_DECADE")
      {
        add_n_per_decade_outputs(n_outputs, n_increments, output_list);
      }
    else
      {
        add_list_outputs(dealii::Utilities::string_to_int(
                           dealii::Utilities::split_string_list(
                             parameter_handler.get("list"))),
                         output_list);
      }
  }
  parameter_handler.leave_subsection();
}

void
FieldOutputParameters::validate(const std::vector<FieldAttributes> &field_attributes,
                                const std::vector<SolveBlock>      &solve_blocks) const
{
  // TODO: Do this later
}

bool
FieldOutputParameters::should_output(unsigned int increment) const
{
  return output_list.contains(increment);
}

void
RestartOutputParameters::predeclare(dealii::ParameterHandler &parameter_handler) const
{
  AssertThrow(false, dealii::ExcNotImplemented());
}

void
RestartOutputParameters::preassign(dealii::ParameterHandler &parameter_handler)
{
  AssertThrow(false, dealii::ExcNotImplemented());
}

void
RestartOutputParameters::declare(dealii::ParameterHandler &parameter_handler,
                                 unsigned int              max_criteria) const
{
  parameter_handler.enter_subsection("checkpoint");
  {
    parameter_handler.declare_entry(
      "load from checkpoint",
      "false",
      dealii::Patterns::Bool(),
      "Whether to load from a checkpoint created during a previous simulation.");

    parameter_handler.declare_entry("directory",
                                    "outputs",
                                    dealii::Patterns::Anything(),
                                    "The name of the output directory.");
    parameter_handler.declare_alias("directory", "folder name");

    parameter_handler.declare_entry("file name",
                                    "solution",
                                    dealii::Patterns::Anything(),
                                    "The prefix of the output files, before the "
                                    "time step and processor info are added.");

    parameter_handler.declare_entry(
      "condition",
      "EQUAL_SPACING",
      dealii::Patterns::Selection("EQUAL_SPACING|LOG_SPACING|N_PER_DECADE|LIST"),
      "The spacing type for outputting the solution fields.");
    parameter_handler.declare_entry(
      "list",
      "0",
      dealii::Patterns::List(dealii::Patterns::Integer(0, INT_MAX), 0, INT_MAX, ","),
      "The list of time steps to output. Used for the LIST type only and must be comma "
      "delimited.");
    parameter_handler.declare_entry("number",
                                    "10",
                                    dealii::Patterns::Integer(0, INT_MAX),
                                    "The number of outputs (or number of outputs "
                                    "per decade for the N_PER_DECADE type).");
  }
  parameter_handler.leave_subsection();
}

void
RestartOutputParameters::assign(dealii::ParameterHandler &parameter_handler,
                                unsigned int              max_criteria)
{
  AssertThrow(false, dealii::ExcNotImplemented());
}

void
RestartOutputParameters::assign(dealii::ParameterHandler &parameter_handler,
                                unsigned int              n_increments,
                                unsigned int              max_criteria)
{
  parameter_handler.enter_subsection("checkpoint");
  {
    load_from_checkpoint = parameter_handler.get_bool("load from checkpoint");
    folder               = parameter_handler.get("directory");
    file_name            = parameter_handler.get("file name");

    std::string  condition = parameter_handler.get("condition");
    unsigned int n_outputs = (unsigned int) (parameter_handler.get_integer("number"));

    if (condition == "EQUAL_SPACING")
      {
        add_equal_spacing_outputs(n_outputs, n_increments, output_list);
      }
    else if (condition == "LOG_SPACING")
      {
        add_log_spacing_outputs(n_outputs, n_increments, output_list);
      }
    else if (condition == "N_PER_DECADE")
      {
        add_n_per_decade_outputs(n_outputs, n_increments, output_list);
      }
    else
      {
        add_list_outputs(dealii::Utilities::string_to_int(
                           dealii::Utilities::split_string_list(
                             parameter_handler.get("list"))),
                         output_list);
      }
  }
  parameter_handler.leave_subsection();
}

void
RestartOutputParameters::validate(const std::vector<FieldAttributes> &field_attributes,
                                  const std::vector<SolveBlock>      &solve_blocks) const
{
  // TODO: Do this later
}

bool
RestartOutputParameters::should_output(unsigned int increment) const
{
  return output_list.contains(increment);
}

void
FieldInputParameters::predeclare(dealii::ParameterHandler &parameter_handler) const
{
  AssertThrow(false, dealii::ExcNotImplemented());
}

void
FieldInputParameters::preassign(dealii::ParameterHandler &parameter_handler)
{
  AssertThrow(false, dealii::ExcNotImplemented());
}

void
FieldInputParameters::declare(dealii::ParameterHandler &parameter_handler,
                              unsigned int              max_criteria) const
{
  for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
    {
      std::string subsection_text = "input file: " + std::to_string(criterion_id);
      parameter_handler.enter_subsection(subsection_text);
      {
        parameter_handler.declare_entry("file name",
                                        "",
                                        dealii::Patterns::Anything(),
                                        "The file name to load from for each variable.");
        parameter_handler.declare_entry("format",
                                        "vtu",
                                        dealii::Patterns::Selection(
                                          "vtk|vtu|vti|pvtu|binary"),
                                        "The type of grid in the file.");
        parameter_handler.declare_entry(
          "file variables",
          "",
          dealii::Patterns::List(dealii::Patterns::Anything(), 0, INT_MAX, ","),
          "The names of the fields in the file.");
        parameter_handler.declare_entry(
          "simulation variables",
          "",
          dealii::Patterns::List(dealii::Patterns::Anything(), 0, INT_MAX, ","),
          "The correspond names of the fields in the simulation.");

        for (const auto &dir : axis_labels)
          {
            const std::string axis {dir};

            parameter_handler.declare_entry(axis + " data points",
                                            "0",
                                            dealii::Patterns::Integer(0, INT_MAX),
                                            "The number of data points in the " + axis +
                                              "-direction.");
          }
      }
      parameter_handler.leave_subsection();
    }

  parameter_handler.declare_entry("read initial conditions from file",
                                  "false",
                                  dealii::Patterns::Bool(),
                                  "Whether to read any initial conditions from file.");
}

void
FieldInputParameters::assign(dealii::ParameterHandler &parameter_handler,
                             unsigned int              max_criteria)
{
  for (unsigned int criterion_id = 0; criterion_id < max_criteria; criterion_id++)
    {
      std::string subsection_text = "input file: " + std::to_string(criterion_id);
      parameter_handler.enter_subsection(subsection_text);
      {
        // Check if the file is specified
        if (!parameter_handler.get("file name").empty())
          {
            InitialConditionFile ic_file;
            ic_file.file_name           = parameter_handler.get("file name");
            ic_file.file_variable_names = dealii::Utilities::split_string_list(
              parameter_handler.get("file variables"));
            ic_file.simulation_variable_names = dealii::Utilities::split_string_list(
              parameter_handler.get("simulation variables"));

            const std::string format = parameter_handler.get("format");
            // TODO: Make this a map like the other enums
            if (format == "vtu")
              {
                ic_file.format =
                  InitialConditionFile::DataFormatType::VTKXMLUnstructuredGrid;
              }
            else if (format == "vtk")
              {
                ic_file.format =
                  InitialConditionFile::DataFormatType::VTKUnstructuredGrid;
              }
            else if (format == "vti")
              {
                ic_file.format = InitialConditionFile::DataFormatType::VTKXMLImageData;
              }
            else if (format == "pvtu")
              {
                ic_file.format =
                  InitialConditionFile::DataFormatType::VTKPXMLUnstructuredGrid;
              }
            else
              {
                ic_file.format = InitialConditionFile::DataFormatType::FlatBinary;
              }

            for (unsigned int k = 0; k < 3; ++k)
              {
                const std::string axis {axis_labels.at(k)};

                ic_file.n_data_points.at(k) = static_cast<unsigned int>(
                  parameter_handler.get_integer(axis + " data points"));
              }

            initial_condition_files.push_back(ic_file);
          }
      }
      parameter_handler.leave_subsection();
    }

  load_from_file = parameter_handler.get_bool("read initial conditions from file");
}

void
FieldInputParameters::validate(const std::vector<FieldAttributes> &field_attributes,
                               const std::vector<SolveBlock>      &solve_blocks) const
{
  // TODO: Do this later
}

PRISMS_PF_END_NAMESPACE
