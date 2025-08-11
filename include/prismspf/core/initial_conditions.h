// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/pde_operator.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>
#include <prismspf/field_input/read_vtk.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

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
    const FieldType                                               &field_type,
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

  std::shared_ptr<const PDEOperator<dim, degree, number>> pde_operator;
};

/**
 * @brief Function for read-in intial conditions. These are only ever used for
 * for explicit time dependent fields and implicit time dependent, as all others are
 * calculated at runtime.
 */
template <unsigned int dim, typename number>
class ReadInitialCondition : public dealii::Function<dim, number>
{
public:
  /**
   * @brief Constructor.
   */
  ReadInitialCondition(const std::string &file_name,
                       std::string        _field_name,
                       const FieldType   &_field_type);

  // NOLINTBEGIN(readability-identifier-length)

  /**
   * @brief Scalar/Vector value.
   */
  void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<number> &value) const override;

  // NOLINTEND(readability-identifier-length)

private:
  std::string field_name;

  FieldType field_type;

  std::shared_ptr<ReadUnstructuredVTK<dim, number>> reader;
};

PRISMS_PF_END_NAMESPACE
