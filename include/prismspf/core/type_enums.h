// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/config.h>

#include <cstdint>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Type of PDE that is being solved.
 */
enum PDEType : std::uint8_t
{
  Constant,
  Explicit,
  Linear,
  Newton
};

/**
 * @brief Tensor rank of the field.
 *
 * Currently, only scalar and vectors are supported.
 */
enum class TensorRank : unsigned int
{
  Undefined = static_cast<unsigned int>(-1),
  Scalar    = 0,
  Vector    = 1
};

/**
 * @brief Symmetry of elastic tensor.
 */
enum ElasticityModel : std::uint8_t
{
  Isotropic,
  Transverse,
  Orthotropic,
  Anisotropic
};

/**
 * @brief Internal classification for types of variable dependencies.
 */
enum DependencyType : int
{
  Normal = -1,
  Change,
  OldOne,
  OldTwo,
  OldThree,
  OldFour
};

/**
 * @brief Solver tolerance type.
 */
enum SolverToleranceType : std::uint8_t
{
  AbsoluteResidual,
  RelativeResidualChange
};

/**
 * @brief Preconditioner type.
 */
enum PreconditionerType : std::uint8_t
{
  None,
  GMG
};

/**
 * @brief Data formats for input initial conditions.
 * LastEntry is used for loop bounds.
 */
enum DataFormatType : std::uint8_t
{
  VTKUnstructuredGrid,
  FlatBinary,
  LastEntry
};

/**
 * @brief Enum to string for ElasticityModel
 */
inline std::string
to_string(ElasticityModel type)
{
  switch (type)
    {
      case ElasticityModel::Isotropic:
        return "Isotropic";
      case ElasticityModel::Transverse:
        return "Transverse";
      case ElasticityModel::Orthotropic:
        return "Orthotropic";
      case ElasticityModel::Anisotropic:
        return "Anisotropic";
      default:
        return "UNKNOWN";
    }
}

/**
 * @brief Enum to string for DependencyType
 */
inline std::string
to_string(DependencyType type)
{
  switch (type)
    {
      case DependencyType::Normal:
        return "Normal";
      case DependencyType::Change:
        return "Change";
      case DependencyType::OldOne:
        return "OldOne";
      case DependencyType::OldTwo:
        return "OldTwo";
      case DependencyType::OldThree:
        return "OldThree";
      case DependencyType::OldFour:
        return "OldFour";
      default:
        return "UNKNOWN";
    }
}

/**
 * @brief Enum to string for SolverToleranceType
 */
inline std::string
to_string(SolverToleranceType type)
{
  switch (type)
    {
      case SolverToleranceType::AbsoluteResidual:
        return "AbsoluteResidual";
      case SolverToleranceType::RelativeResidualChange:
        return "RelativeResidualChange";
      default:
        return "UNKNOWN";
    }
}

/**
 * @brief Enum to string for PreconditionerType
 */
inline std::string
to_string(PreconditionerType type)
{
  switch (type)
    {
      case PreconditionerType::None:
        return "None";
      case PreconditionerType::GMG:
        return "GMG";
      default:
        return "UNKNOWN";
    }
}

/**
 * @brief Enum to string for DataFormatType
 */
inline std::string
to_string(DataFormatType type)
{
  switch (type)
    {
      case DataFormatType::VTKUnstructuredGrid:
        return "vtk_unstructured_grid";
      case DataFormatType::FlatBinary:
        return "flat_binary";
      default:
        return "unknown";
    }
}

PRISMS_PF_END_NAMESPACE
