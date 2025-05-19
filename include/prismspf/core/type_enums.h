// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/config.h>

#include <cstdint>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Type of field. Currently, the only support fields are scalar and vector.
 */
enum FieldType : std::uint8_t
{
  UNDEFINED_FIELD,
  SCALAR,
  VECTOR
};

/**
 * \brief Type of PDE that is being solved.
 */
enum PDEType : std::uint8_t
{
  UNDEFINED_PDE,
  EXPLICIT_TIME_DEPENDENT,
  IMPLICIT_TIME_DEPENDENT,
  TIME_INDEPENDENT,
  AUXILIARY,
  CONSTANT
};

/**
 * \brief Type of solve.
 */
enum SolveType : std::uint8_t
{
  EXPLICIT_rhs,
  NONEXPLICIT_rhs,
  NONEXPLICIT_lhs,
  POSTPROCESS
};

/**
 * \brief Symmetry of elastic tensor.
 */
enum ElasticityModel : std::uint8_t
{
  ISOTROPIC,
  TRANSVERSE,
  ORTHOTROPIC,
  ANISOTROPIC
};

/**
 * \brief Internal classification of combined field and solve types. There are six
 * different types of solve that are possible. For EXPLICIT solves, all fields of that
 * type can be solved concurrently. For NONEXPLICIT_LINEAR, NONEXPLICIT_SELF_NONLINEAR,
 * and NONEXPLICIT_AUXILIARY, these must be solved sequentially and wrapped in
 * conditionals in the user implmentation. For NONEXPLICIT_CO_NONLINEAR, there are at
 * least 2 fields that are nonlinear together, as opposed to NONEXPLICIT_SELF_NONLINEAR,
 * which must be solved at the same time. A simply case for this is the steady-state
 * Cahn-Hilliard equation. Finally, for EXPLICIT_POSTPROCESS and EXPLICIT_CONSTANT, they
 * are more or less the same as EXPLICIT.
 */
enum FieldSolveType : std::uint8_t
{
  UNDEFINED_SOLVE,
  EXPLICIT_CONSTANT,
  EXPLICIT,
  NONEXPLICIT_LINEAR,
  NONEXPLICIT_SELF_NONLINEAR,
  NONEXPLICIT_AUXILIARY,
  NONEXPLICIT_CO_NONLINEAR,
  EXPLICIT_POSTPROCESS,
};

/**
 * \brief Internal classification for types of variable dependencies.
 */
enum DependencyType : std::uint8_t
{
  NORMAL,
  CHANGE,
  OLD_1,
  OLD_2,
  OLD_3,
  OLD_4
};

/**
 * \brief Solver tolerance type.
 */
enum SolverToleranceType : std::uint8_t
{
  ABSOLUTE_RESIDUAL,
  RELATIVE_RESIDUAL_CHANGE
};

/**
 * \brief Preconditioner type.
 */
enum PreconditionerType : std::uint8_t
{
  NONE,
  GMG
};

/**
 * \brief Enum to string for FieldType
 */
inline std::string
to_string(FieldType type)

{
  switch (type)
    {
      case FieldType::UNDEFINED_FIELD:
        return "UNDEFINED_FIELD";
      case FieldType::SCALAR:
        return "SCALAR_FIELD";
      case FieldType::VECTOR:
        return "VECTOR_FIELD";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for PDEType
 */
inline std::string
to_string(PDEType type)
{
  switch (type)
    {
      case PDEType::UNDEFINED_PDE:
        return "UNDEFINED_PDE";
      case PDEType::EXPLICIT_TIME_DEPENDENT:
        return "EXPLICIT_TIME_DEPENDENT";
      case PDEType::IMPLICIT_TIME_DEPENDENT:
        return "IMPLICIT_TIME_DEPENDENT";
      case PDEType::TIME_INDEPENDENT:
        return "TIME_INDEPENDENT";
      case PDEType::AUXILIARY:
        return "AUXILIARY";
      case PDEType::CONSTANT:
        return "CONSTANT";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for SolveType
 */
inline std::string
to_string(SolveType type)
{
  switch (type)
    {
      case SolveType::EXPLICIT_rhs:
        return "EXPLICIT_rhs";
      case SolveType::NONEXPLICIT_rhs:
        return "NONEXPLICIT_rhs";
      case SolveType::NONEXPLICIT_lhs:
        return "NONEXPLICIT_lhs";
      case SolveType::POSTPROCESS:
        return "POSTPROCESS";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for ElasticityModel
 */
inline std::string
to_string(ElasticityModel type)
{
  switch (type)
    {
      case ElasticityModel::ISOTROPIC:
        return "ISOTROPIC";
      case ElasticityModel::TRANSVERSE:
        return "TRANSVERSE";
      case ElasticityModel::ORTHOTROPIC:
        return "ORTHOTROPIC";
      case ElasticityModel::ANISOTROPIC:
        return "ANISOTROPIC";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for FieldSolveType
 */
inline std::string
to_string(FieldSolveType type)
{
  switch (type)
    {
      case FieldSolveType::UNDEFINED_SOLVE:
        return "UNDEFINED_SOLVE";
      case FieldSolveType::EXPLICIT:
        return "EXPLICIT";
      case FieldSolveType::NONEXPLICIT_LINEAR:
        return "NONEXPLICIT_LINEAR";
      case FieldSolveType::NONEXPLICIT_SELF_NONLINEAR:
        return "NONEXPLICIT_SELF_NONLINEAR";
      case FieldSolveType::NONEXPLICIT_AUXILIARY:
        return "NONEXPLICIT_AUXILIARY";
      case FieldSolveType::NONEXPLICIT_CO_NONLINEAR:
        return "NONEXPLICIT_CO_NONLINEAR";
      case FieldSolveType::EXPLICIT_POSTPROCESS:
        return "EXPLICIT_POSTPROCESS";
      case FieldSolveType::EXPLICIT_CONSTANT:
        return "EXPLICIT_CONSTANT";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for DependencyType
 */
inline std::string
to_string(DependencyType type)
{
  switch (type)
    {
      case DependencyType::NORMAL:
        return "NORMAL";
      case DependencyType::CHANGE:
        return "CHANGE";
      case DependencyType::OLD_1:
        return "OLD_1";
      case DependencyType::OLD_2:
        return "OLD_2";
      case DependencyType::OLD_3:
        return "OLD_3";
      case DependencyType::OLD_4:
        return "OLD_4";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for SolverToleranceType
 */
inline std::string
to_string(SolverToleranceType type)
{
  switch (type)
    {
      case SolverToleranceType::ABSOLUTE_RESIDUAL:
        return "ABSOLUTE_RESIDUAL";
      case SolverToleranceType::RELATIVE_RESIDUAL_CHANGE:
        return "RELATIVE_RESIDUAL_CHANGE";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for PreconditionerType
 */
inline std::string
to_string(PreconditionerType type)
{
  switch (type)
    {
      case PreconditionerType::NONE:
        return "NONE";
      case PreconditionerType::GMG:
        return "GMG";
      default:
        return "UNKNOWN";
    }
}

PRISMS_PF_END_NAMESPACE