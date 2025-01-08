#ifndef type_enums_h
#define type_enums_h

#include <cstdint>
#include <string>

/**
 * \brief Type of field. Currently, the only support fields are scalar and vector.
 */
enum fieldType : std::uint8_t
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
enum solveType : std::uint8_t
{
  EXPLICIT_RHS,
  NONEXPLICIT_RHS,
  NONEXPLICIT_LHS,
  POSTPROCESS
};

/**
 * \brief Symmetry of elastic tensor.
 */
enum elasticityModel : std::uint8_t
{
  ISOTROPIC,
  TRANSVERSE,
  ORTHOTROPIC,
  ANISOTROPIC
};

/**
 * \brief Type of boundary condition
 */
enum boundaryType : std::uint8_t
{
  UNDEFINED_BOUNDARY,
  NATURAL,
  DIRICHLET,
  PERIODIC,
  NEUMANN,
  NON_UNIFORM_DIRICHLET,
  NON_UNIFORM_NEUMANN
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
enum fieldSolveType : std::uint8_t
{
  UNDEFINED_SOLVE,
  EXPLICIT,
  NONEXPLICIT_LINEAR,
  NONEXPLICIT_SELF_NONLINEAR,
  NONEXPLICIT_AUXILIARY,
  NONEXPLICIT_CO_NONLINEAR,
  EXPLICIT_POSTPROCESS,
  EXPLICIT_CONSTANT
};

/**
 * \brief Internal classification for types of variable dependencies.
 */
enum dependencyType : std::uint8_t
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
enum solverToleranceType : std::uint8_t
{
  ABSOLUTE_RESIDUAL,
  RELATIVE_RESIDUAL_CHANGE
};

/**
 * \brief Preconditioner type.
 */
enum preconditionerType : std::uint8_t
{
  NONE,
  GMG
};

/**
 * \brief Enum to string for fieldType
 */
inline std::string
to_string(fieldType type)

{
  switch (type)
    {
      case fieldType::UNDEFINED_FIELD:
        return "UNDEFINED_FIELD";
      case fieldType::SCALAR:
        return "SCALAR_FIELD";
      case fieldType::VECTOR:
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
 * \brief Enum to string for solveType
 */
inline std::string
to_string(solveType type)
{
  switch (type)
    {
      case solveType::EXPLICIT_RHS:
        return "EXPLICIT_RHS";
      case solveType::NONEXPLICIT_RHS:
        return "NONEXPLICIT_RHS";
      case solveType::NONEXPLICIT_LHS:
        return "NONEXPLICIT_LHS";
      case solveType::POSTPROCESS:
        return "POSTPROCESS";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for elasticityModel
 */
inline std::string
to_string(elasticityModel type)
{
  switch (type)
    {
      case elasticityModel::ISOTROPIC:
        return "ISOTROPIC";
      case elasticityModel::TRANSVERSE:
        return "TRANSVERSE";
      case elasticityModel::ORTHOTROPIC:
        return "ORTHOTROPIC";
      case elasticityModel::ANISOTROPIC:
        return "ANISOTROPIC";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for boundaryType
 */
inline std::string
to_string(boundaryType type)
{
  switch (type)
    {
      case boundaryType::UNDEFINED_BOUNDARY:
        return "UNDEFINED_BOUNDARY";
      case boundaryType::NATURAL:
        return "NATURAL";
      case boundaryType::DIRICHLET:
        return "DIRICHLET";
      case boundaryType::PERIODIC:
        return "PERIODIC";
      case boundaryType::NEUMANN:
        return "NEUMANN";
      case boundaryType::NON_UNIFORM_DIRICHLET:
        return "NON_UNIFORM_DIRICHLET";
      case boundaryType::NON_UNIFORM_NEUMANN:
        return "NON_UNIFORM_NEUMANN";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for fieldSolveType
 */
inline std::string
to_string(fieldSolveType type)
{
  switch (type)
    {
      case fieldSolveType::UNDEFINED_SOLVE:
        return "UNDEFINED_SOLVE";
      case fieldSolveType::EXPLICIT:
        return "EXPLICIT";
      case fieldSolveType::NONEXPLICIT_LINEAR:
        return "NONEXPLICIT_LINEAR";
      case fieldSolveType::NONEXPLICIT_SELF_NONLINEAR:
        return "NONEXPLICIT_SELF_NONLINEAR";
      case fieldSolveType::NONEXPLICIT_AUXILIARY:
        return "NONEXPLICIT_AUXILIARY";
      case fieldSolveType::NONEXPLICIT_CO_NONLINEAR:
        return "NONEXPLICIT_CO_NONLINEAR";
      case fieldSolveType::EXPLICIT_POSTPROCESS:
        return "EXPLICIT_POSTPROCESS";
      case fieldSolveType::EXPLICIT_CONSTANT:
        return "EXPLICIT_CONSTANT";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for dependencyType
 */
inline std::string
to_string(dependencyType type)
{
  switch (type)
    {
      case dependencyType::NORMAL:
        return "NORMAL";
      case dependencyType::CHANGE:
        return "CHANGE";
      case dependencyType::OLD_1:
        return "OLD_1";
      case dependencyType::OLD_2:
        return "OLD_2";
      case dependencyType::OLD_3:
        return "OLD_3";
      case dependencyType::OLD_4:
        return "OLD_4";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for solverToleranceType
 */
inline std::string
to_string(solverToleranceType type)
{
  switch (type)
    {
      case solverToleranceType::ABSOLUTE_RESIDUAL:
        return "ABSOLUTE_RESIDUAL";
      case solverToleranceType::RELATIVE_RESIDUAL_CHANGE:
        return "RELATIVE_RESIDUAL_CHANGE";
      default:
        return "UNKNOWN";
    }
}

/**
 * \brief Enum to string for preconditionerType
 */
inline std::string
to_string(preconditionerType type)
{
  switch (type)
    {
      case preconditionerType::NONE:
        return "NONE";
      case preconditionerType::GMG:
        return "GMG";
      default:
        return "UNKNOWN";
    }
}

#endif
