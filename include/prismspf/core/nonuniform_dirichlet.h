// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Forward declaration of user-facing implementation
 */
template <int dim>
class customNonuniformDirichlet;

/**
 * \brief Function for user-implemented nonuniform dirichlet boundary condition.
 */
template <int dim, fieldType field_type = fieldType::SCALAR>
class nonuniformDirichlet : public dealii::Function<dim, double>
{
public:
  /**
   * \brief Constructor.
   */
  nonuniformDirichlet(const unsigned int             &_index,
                      const unsigned int             &_boundary_id,
                      const userInputParameters<dim> &_user_inputs);

  // NOLINTBEGIN(readability-identifier-length, readability-avoid-const-params-in-decls)

  /**
   * \brief Scalar value.
   */
  double
  value(const dealii::Point<dim> &p, const unsigned int component = 0) const override;

  /**
   * \brief Vector value.
   */
  void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<double> &value) const override;

  // NOLINTEND(readability-identifier-length, readability-avoid-const-params-in-decls)

private:
  unsigned int index;

  unsigned int boundary_id;

  const userInputParameters<dim> *user_inputs;

  customNonuniformDirichlet<dim> custom_nonuniform_dirichlet;
};

template <int dim, fieldType field_type>
nonuniformDirichlet<dim, field_type>::nonuniformDirichlet(
  const unsigned int             &_index,
  const unsigned int             &_boundary_id,
  const userInputParameters<dim> &_user_inputs)
  : dealii::Function<dim>((field_type == fieldType::VECTOR) ? dim : 1)
  , index(_index)
  , boundary_id(_boundary_id)
  , user_inputs(&_user_inputs)
{}

// NOLINTBEGIN(readability-identifier-length)

template <int dim, fieldType field_type>
inline double
nonuniformDirichlet<dim, field_type>::value(
  const dealii::Point<dim>           &p,
  [[maybe_unused]] const unsigned int component) const
{
  // Initialize passed variables to zero
  double                 temp_scalar_value = 0.0;
  dealii::Vector<double> temp_vector_value(dim);

  // Pass variables to user-facing function to evaluate
  custom_nonuniform_dirichlet.set_nonuniform_dirichlet(index,
                                                       boundary_id,
                                                       0,
                                                       p,
                                                       temp_scalar_value,
                                                       temp_vector_value(0),
                                                       *user_inputs);

  return temp_scalar_value;
}

template <int dim, fieldType field_type>
inline void
nonuniformDirichlet<dim, field_type>::vector_value(const dealii::Point<dim> &p,
                                                   dealii::Vector<double>   &value) const
{
  // Initialize passed variables to zero
  double                 temp_scalar_value = 0.0;
  dealii::Vector<double> temp_vector_value(dim);

  // Pass variables to user-facing function to evaluate
  for (unsigned int i = 0; i < dim; i++)
    {
      custom_nonuniform_dirichlet.set_nonuniform_dirichlet(index,
                                                           boundary_id,
                                                           i,
                                                           p,
                                                           temp_scalar_value,
                                                           temp_vector_value(i),
                                                           *user_inputs);
    }

  value = temp_vector_value;
}

// NOLINTEND(readability-identifier-length)

/**
 * \brief User-facing implementation of nonuniform boundary conditions
 */
template <int dim>
class customNonuniformDirichlet
{
public:
  /**
   * \brief Constructor.
   */
  customNonuniformDirichlet() = default;

  /**
   * \brief Function that passes the value/vector and point that are set in the nonuniform
   * dirichlet.
   */
  void
  set_nonuniform_dirichlet(const unsigned int             &index,
                           const unsigned int             &boundary_id,
                           const unsigned int             &component,
                           const dealii::Point<dim>       &point,
                           double                         &scalar_value,
                           double                         &vector_component_value,
                           const userInputParameters<dim> &user_inputs) const;
};

PRISMS_PF_END_NAMESPACE
