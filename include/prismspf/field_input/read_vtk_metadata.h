// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/config.h>

#include <optional>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

#ifdef PRISMS_PF_WITH_VTK

/**
 * @brief Simulated time and increment metadata read from a VTK file's FieldData.
 */
struct VTKTimeCycleMetadata
{
  double       time      = 0.0;
  unsigned int increment = 0;
};

/**
 * @brief Attempt to read TIME/CYCLE FieldData metadata from a legacy VTK file.
 *
 * Returns std::nullopt if the file lacks TIME/CYCLE FieldData (not an error).
 *
 * Only supported for the legacy VTK ("vtk") unstructured grid format.
 */
std::optional<VTKTimeCycleMetadata>
read_vtk_time_cycle_metadata(const std::string &file_name);

#endif

PRISMS_PF_END_NAMESPACE
