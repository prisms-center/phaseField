// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <prismspf/core/types.h>

#include <prismspf/utilities/utilities.h>
#include <prismspf/user_inputs/spatial_discretization.h>
#include <prismspf/user_inputs/load_initial_condition_parameters.h>

#include <prismspf/field_input/read_field_base.h>
#include <prismspf/field_input/read_vtk.h>


PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Factory function to create appropriate reader based on input file type
 * not a member of ReadFieldBase to avoid redundant template instantiations
 */

enum class Type { ReadUnstructuredVTK };

template <unsigned int dim, typename number>
std::shared_ptr<ReadFieldBase<dim, number>>
create_reader(const InitialConditionFile &ic_file,
              const SpatialDiscretization<dim> &spatial_discretization)
{
  const static std::unordered_map<std::string,int> string_to_case
    {
      {"vtk_unstructured_grid",1},
    };
  switch (string_to_case.contains(ic_file.dataset_format) ? string_to_case.at(ic_file.dataset_format) : 0)
   {
      case 1:
          return std::make_shared<ReadUnstructuredVTK<dim, number>>(ic_file.filename);
      default:
          AssertThrow(false, dealii::ExcMessage("Unsupported dataset format: " + ic_file.dataset_format));
    }
}

PRISMS_PF_END_NAMESPACE