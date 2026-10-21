// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/user_inputs/io_parameters.h>

#include <prismspf/config.h>

#include "prismspf/user_inputs/parameter_base.h"

PRISMS_PF_BEGIN_NAMESPACE

void
FieldOutputParameters::declare(dealii::ParameterHandler &parameter_handler,
                               unsigned int              n_subsections)
{
  parameter_handler.enter_subsection("output");
  {
    declare_entry(parameter_handler,
                  "file type",
                  "vtu",
                  Patterns::Selection("vtu|vtk|pvtu|xdmf"),
                  "The output file type.",
                  {"type", "file extension"});
    declare_entry(parameter_handler,
                  "subdivisions",
                  "0",
                  Patterns::UnsignedInteger(),
                  "The number of subdivisions to apply to the mesh when building output "
                  "patches. If 0, the degree is used.",
                  {"n patches", "number of patches", "num patches"});
    declare_entry(
      parameter_handler,
      "compression level",
      "default",
      Patterns::Selection("default|speed|size|none"),
      "The compression level for output. `speed` is fast but results in larger files, "
      "`size` is slow but gives the best compression, and `none` skips compression "
      "entirely.",
      {"compression"});
    declare_entry(parameter_handler,
                  "directory",
                  "solutions",
                  Patterns::Anything(),
                  "The name of the output directory.",
                  {"folder name"});
    declare_entry(parameter_handler,
                  "file name",
                  "solution",
                  Patterns::Anything(),
                  "The prefix of the output files, before the time step and processor "
                  "info are added.",
                  {"base file name", "base name"});
    declare_entry(parameter_handler,
                  "condition",
                  "EQUAL_SPACING",
                  Patterns::Selection("EQUAL_SPACING|LOG_SPACING|N_PER_DECADE"),
                  "The spacing type for outputting the solution fields.");
    declare_entry(parameter_handler,
                  "number",
                  "10",
                  Patterns::UnsignedInteger(),
                  "The number of outputs for EQUAL_SPACING and for LOG_SPACING or number "
                  "of outputs per decade for N_PER_DECADE.");
    declare_entry(parameter_handler,
                  "list",
                  "0",
                  Patterns::UnsignedIntegerList(),
                  "A list of increments to output on. Note that this is in addition to "
                  "the `condition` and `number` above.");
    declare_entry(
      parameter_handler,
      "variables",
      "",
      Patterns::StringList(),
      "The list of the fields to output. Additionally, for the output of left-hand side "
      "and old fields, they must follow the same delimiters that are used in dependency "
      "sets.");
  }
  parameter_handler.leave_subsection();
}

void
FieldOutputParameters::assign(dealii::ParameterHandler &parameter_handler,
                              unsigned int              n_subsections)
{
  AssertThrow(false, dealii::ExcNotImplemented());
}

void
FieldOutputParameters::assign(dealii::ParameterHandler &parameter_handler,
                              unsigned int              n_increments,
                              unsigned int              n_subsections)
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
    add_list_outputs(dealii::Utilities::string_to_int(
                       dealii::Utilities::split_string_list(
                         parameter_handler.get("list"))),
                     output_list);
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
RestartOutputParameters::declare(dealii::ParameterHandler &parameter_handler,
                                 unsigned int              n_subsections)
{
  parameter_handler.enter_subsection("checkpoint");
  {
    declare_entry(
      parameter_handler,
      "load from checkpoint",
      "false",
      Patterns::Bool(),
      "Whether to load from a checkpoint created during a previous simulation.");
    declare_entry(parameter_handler,
                  "directory",
                  "solutions",
                  Patterns::Anything(),
                  "The name of the output directory.",
                  {"folder name"});
    declare_entry(parameter_handler,
                  "file name",
                  "solution",
                  Patterns::Anything(),
                  "The prefix of the output files, before the time step and processor "
                  "info are added.");
    declare_entry(parameter_handler,
                  "condition",
                  "EQUAL_SPACING",
                  Patterns::Selection("EQUAL_SPACING|LOG_SPACING|N_PER_DECADE"),
                  "The spacing type for outputting the solution fields.");
    declare_entry(parameter_handler,
                  "list",
                  "0",
                  Patterns::UnsignedIntegerList(),
                  "The list of time steps to output.");
    declare_entry(parameter_handler,
                  "number",
                  "10",
                  Patterns::UnsignedInteger(),
                  "The number of outputs (or number of outputs "
                  "per decade for the N_PER_DECADE type).");
  }
  parameter_handler.leave_subsection();
}

void
RestartOutputParameters::assign(dealii::ParameterHandler &parameter_handler,
                                unsigned int              n_subsections)
{
  AssertThrow(false, dealii::ExcNotImplemented());
}

void
RestartOutputParameters::assign(dealii::ParameterHandler &parameter_handler,
                                unsigned int              n_increments,
                                unsigned int              n_subsections)
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
FieldInputParameters::declare(dealii::ParameterHandler &parameter_handler,
                              unsigned int              n_subsections)
{
  for (unsigned int criterion_id = 0; criterion_id < n_subsections; criterion_id++)
    {
      std::string subsection_text = "input file: " + std::to_string(criterion_id);
      parameter_handler.enter_subsection(subsection_text);
      {
        declare_entry(parameter_handler,
                      "file name",
                      "",
                      dealii::Patterns::Anything(),
                      "The file name to load from for each variable.");
        declare_entry(parameter_handler,
                      "format",
                      "vtu",
                      dealii::Patterns::Selection("vtk|vtu|vti|pvtu|binary"),
                      "The type of grid in the file.");
        declare_entry(parameter_handler,
                      "file variables",
                      "",
                      Patterns::StringList(),
                      "The names of the fields in the file.");
        declare_entry(parameter_handler,
                      "simulation variables",
                      "",
                      Patterns::StringList(),
                      "The correspond names of the fields in the simulation.");

        for (const auto &dir : axis_labels)
          {
            const std::string axis {dir};

            declare_entry(parameter_handler,
                          axis + " data points",
                          "0",
                          Patterns::UnsignedInteger(),
                          "The number of data points in the " + axis + "-direction.");
          }
      }
      parameter_handler.leave_subsection();
    }
  declare_entry(parameter_handler,
                "read initial conditions from file",
                "false",
                Patterns::Bool(),
                "Whether to read any initial conditions from file.");
}

void
FieldInputParameters::assign(dealii::ParameterHandler &parameter_handler,
                             unsigned int              n_subsections)
{
  for (unsigned int criterion_id = 0; criterion_id < n_subsections; criterion_id++)
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
