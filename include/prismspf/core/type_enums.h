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
  ExplicitTimeDependent,
  ImplicitTimeDependent,
  TimeIndependent,
  Auxiliary,
  Constant
};

/**
 * @brief Type of solve.
 */
enum SolveType : std::uint8_t
{
  ExplicitRHS,
  NonexplicitRHS,
  NonexplicitLHS,
  Postprocess
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
 * @brief Internal classification of combined field and solve types. There are six
 * different types of solve that are possible. For Explicit solves, all fields of that
 * type can be solved concurrently. For NonexplicitLinear, NonexplicitSelfnonlinear,
 * and NonexplicitAuxiliary, these must be solved sequentially and wrapped in
 * conditionals in the user implmentation. For NonexplicitCononlinear, there are at
 * least 2 fields that are nonlinear together, as opposed to NonexplicitSelfnonlinear,
 * which must be solved at the same time. A simply case for this is the steady-state
 * Cahn-Hilliard equation. Finally, for ExplicitPostprocess and ExplicitConstant, they
 * are more or less the same as Explicit.
 */
enum FieldSolveType : std::uint8_t
{
  ExplicitConstant,
  Explicit,
  NonexplicitLinear,
  NonexplicitSelfnonlinear,
  NonexplicitAuxiliary,
  NonexplicitCononlinear,
  ExplicitPostprocess,
};

/**
 * @brief Internal classification for types of variable dependencies.
 */
enum DependencyType : std::uint8_t
{
  Normal,
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
 * @brief Enum to string for PDEType
 */
inline std::string
to_string(PDEType type)
{
  switch (type)
    {
      case PDEType::ExplicitTimeDependent:
        return "ExplicitTimeDependent";
      case PDEType::ImplicitTimeDependent:
        return "ImplicitTimeDependent";
      case PDEType::TimeIndependent:
        return "TimeIndependent";
      case PDEType::Auxiliary:
        return "Auxiliary";
      case PDEType::Constant:
        return "Constant";
      default:
        return "UNKNOWN";
    }
}

/**
 * @brief Enum to string for SolveType
 */
inline std::string
to_string(SolveType type)
{
  switch (type)
    {
      case SolveType::ExplicitRHS:
        return "ExplicitRHS";
      case SolveType::NonexplicitRHS:
        return "NonexplicitRHS";
      case SolveType::NonexplicitLHS:
        return "NonexplicitLHS";
      case SolveType::Postprocess:
        return "Postprocess";
      default:
        return "UNKNOWN";
    }
}

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
 * @brief Enum to string for FieldSolveType
 */
inline std::string
to_string(FieldSolveType type)
{
  switch (type)
    {
      case FieldSolveType::Explicit:
        return "Explicit";
      case FieldSolveType::NonexplicitLinear:
        return "NonexplicitLinear";
      case FieldSolveType::NonexplicitSelfnonlinear:
        return "NonexplicitSelfnonlinear";
      case FieldSolveType::NonexplicitAuxiliary:
        return "NonexplicitAuxiliary";
      case FieldSolveType::NonexplicitCononlinear:
        return "NonexplicitCononlinear";
      case FieldSolveType::ExplicitPostprocess:
        return "ExplicitPostprocess";
      case FieldSolveType::ExplicitConstant:
        return "ExplicitConstant";
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
