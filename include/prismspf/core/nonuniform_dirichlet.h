// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef nonuniform_dirichlet_h
#define nonuniform_dirichlet_h

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/config.h>
#include <prismspf/core/type_enums.h>

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
  nonuniformDirichlet(const unsigned int &_index, const unsigned int &_boundary_id);

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

private:
  const unsigned int index;

  const unsigned int boundary_id;

  customNonuniformDirichlet<dim> custom_nonuniform_dirichlet;
};

template <int dim, fieldType field_type>
nonuniformDirichlet<dim, field_type>::nonuniformDirichlet(
  const unsigned int &_index,
  const unsigned int &_boundary_id)
  : dealii::Function<dim>((field_type == fieldType::VECTOR) ? dim : 1)
  , index(_index)
  , boundary_id(_boundary_id)
{}

template <int dim, fieldType field_type>
inline double
nonuniformDirichlet<dim, field_type>::value(
  const dealii::Point<dim>           &p,
  [[maybe_unused]] const unsigned int component) const
{
  // Initialize passed variables to zero
  double                 scalar_value = 0.0;
  dealii::Vector<double> vector_value(dim);

  // Pass variables to user-facing function to evaluate
  custom_nonuniform_dirichlet
    .set_nonuniform_dirichlet(index, boundary_id, 0, p, scalar_value, vector_value(0));

  return scalar_value;
}

template <int dim, fieldType field_type>
inline void
nonuniformDirichlet<dim, field_type>::vector_value(const dealii::Point<dim> &p,
                                                   dealii::Vector<double>   &value) const
{
  // Initialize passed variables to zero
  double                 scalar_value = 0.0;
  dealii::Vector<double> vector_value(dim);

  // Pass variables to user-facing function to evaluate
  for (unsigned int i = 0; i < dim; i++)
    {
      custom_nonuniform_dirichlet.set_nonuniform_dirichlet(index,
                                                           boundary_id,
                                                           i,
                                                           p,
                                                           scalar_value,
                                                           vector_value(i));
    }

  value = vector_value;
}

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
  set_nonuniform_dirichlet(const unsigned int       &index,
                           const unsigned int       &boundary_id,
                           const unsigned int       &component,
                           const dealii::Point<dim> &point,
                           double                   &scalar_value,
                           double                   &vector_component_value
                           const userInputParameters<dim> &user_inputs) const;
};

PRISMS_PF_END_NAMESPACE

#endif