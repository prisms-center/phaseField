// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/field_input/read_vtk_metadata.h>

#include <prismspf/config.h>

#ifdef PRISMS_PF_WITH_VTK

#  include <filesystem>
#  include <vtkDataArray.h>
#  include <vtkFieldData.h>
#  include <vtkNew.h>
#  include <vtkUnstructuredGrid.h>
#  include <vtkUnstructuredGridReader.h>

PRISMS_PF_BEGIN_NAMESPACE

std::optional<VTKTimeCycleMetadata>
read_vtk_time_cycle_metadata(const std::string &file_name)
{
  if (!std::filesystem::exists(file_name))
    {
      return std::nullopt;
    }

  vtkNew<vtkUnstructuredGridReader> reader;
  reader->SetFileName(file_name.c_str());
  reader->Update();

  vtkUnstructuredGrid *output = reader->GetOutput();
  if (output == nullptr)
    {
      return std::nullopt;
    }

  vtkFieldData *field_data = output->GetFieldData();
  if (field_data == nullptr)
    {
      return std::nullopt;
    }

  vtkDataArray *time_array  = field_data->GetArray("TIME");
  vtkDataArray *cycle_array = field_data->GetArray("CYCLE");
  if (time_array == nullptr || cycle_array == nullptr)
    {
      return std::nullopt;
    }

  VTKTimeCycleMetadata metadata;
  metadata.time      = time_array->GetComponent(0, 0);
  metadata.increment = static_cast<unsigned int>(cycle_array->GetComponent(0, 0));
  return metadata;
}

PRISMS_PF_END_NAMESPACE

#endif
