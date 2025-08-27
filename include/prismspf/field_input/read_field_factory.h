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

#if DEAL_II_VERSION_MAJOR >= 9 && DEAL_II_VERSION_MINOR >= 7
#  include <deal.II/base/exception_macros.h>
#else
#  include <deal.II/base/exceptions.h>
#endif

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
  switch (ic_file.dataset_format)
   {
      case DataFormatType::VTKUnstructuredGrid:
        return std::make_shared<ReadUnstructuredVTK<dim, number>>(ic_file.filename);
      case DataFormatType::FlatBinary:
        return std::make_shared<ReadUnstructuredVTK<dim, number>>(ic_file.filename);
      default:
          AssertThrow(false, UnreachableCode());
    }
}

PRISMS_PF_END_NAMESPACE