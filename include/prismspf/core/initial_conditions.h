// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

struct InitialConditionFile;

template <unsigned int dim>
struct SpatialDiscretization;

template <unsigned int dim, unsigned int degree, typename number>
class PDEOperator;

template <unsigned int dim, typename number>
class ReadFieldBase;

/**
 * @brief Function for user-implemented initial conditions. These are only ever calculated
 * for explicit time dependent fields and implicit time dependent, as all others are
 * calculated at runtime.
 */
template <unsigned int dim, unsigned int degree, typename number>
class InitialCondition : public dealii::Function<dim, number>
{
public:
  /**
   * @brief Constructor.
   */
  InitialCondition(
    const unsigned int                                            &_index,
    const FieldInfo::TensorRank                                   &_field_type,
    const std::shared_ptr<const PDEOperator<dim, degree, number>> &_pde_operator);

  // NOLINTBEGIN(readability-identifier-length)

  /**
   * @brief Scalar/Vector value.
   */
  void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<number> &value) const override;

  // NOLINTEND(readability-identifier-length)

private:
  unsigned int index;

  FieldInfo::TensorRank field_type;

  std::shared_ptr<const PDEOperator<dim, degree, number>> pde_operator;
};

/**
 * @brief Function for read-in of initial conditions.
 */
template <unsigned int dim, typename number>
class ReadInitialCondition : public dealii::Function<dim, number>
{
public:
  /**
   * @brief Constructor.
   */
  ReadInitialCondition(std::string                       _field_name,
                       const FieldInfo::TensorRank      &_field_type,
                       const InitialConditionFile       &ic_file,
                       const SpatialDiscretization<dim> &spatial_discretization);

  // NOLINTBEGIN(readability-identifier-length)

  /**
   * @brief Scalar/Vector value.
   */
  void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<number> &value) const override;

  // NOLINTEND(readability-identifier-length)

private:
  std::string field_name;

  FieldInfo::TensorRank field_type;

  std::shared_ptr<ReadFieldBase<dim, number>> reader;
};

PRISMS_PF_END_NAMESPACE
