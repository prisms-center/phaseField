// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>

#include <memory>
#include <string>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
InitialCondition<dim, degree>::InitialCondition(
  const unsigned int                                            &_index,
  const FieldType                                               &field_type,
  const std::shared_ptr<const PDEOperator<dim, degree, double>> &_pde_operator)
  : dealii::Function<dim>((field_type == FieldType::Vector) ? dim : 1)
  , index(_index)
  , pde_operator(_pde_operator)
{}

// NOLINTBEGIN(readability-identifier-length)

template <unsigned int dim, unsigned int degree>
void
InitialCondition<dim, degree>::vector_value(const dealii::Point<dim> &p,
                                            dealii::Vector<double>   &value) const
{
  // Initialize passed variables to zero
  dealii::Vector<double> vector_value(dim);

  // Pass variables to user-facing function to evaluate
  // TODO (landinjm): This is redoing the scalar field calculation.
  for (unsigned int i = 0; i < dim; i++)
    {
      pde_operator->set_initial_condition(index, i, p, vector_value[0], vector_value[i]);
    }

  value = vector_value;
}

// NOLINTEND(readability-identifier-length)

template <unsigned int dim>
ReadInitialCondition<dim>::ReadInitialCondition(const std::string &file_name,
                                                std::string        _field_name,
                                                const FieldType   &_field_type)
  : dealii::Function<dim>((_field_type == FieldType::Vector) ? dim : 1)
  , field_name(std::move(_field_name))
  , field_type(_field_type)
  , reader(std::make_shared<ReadUnstructuredVTK<dim>>(file_name))
{}

// NOLINTBEGIN(readability-identifier-length)

template <unsigned int dim>
void
ReadInitialCondition<dim>::vector_value(const dealii::Point<dim> &p,
                                        dealii::Vector<double>   &value) const
{
  // Initialize passed variables to zero
  dealii::Vector<double> vector_value(dim);

  // Get the value of the field at the point
  if (field_type == FieldType::Scalar)
    {
      vector_value[0] = reader->get_scalar_value(p, field_name);
    }
  else if (field_type == FieldType::Vector)
    {
      vector_value = reader->get_vector_value(p, field_name);
    }

  value = vector_value;
}

// NOLINTEND(readability-identifier-length)

#include "core/initial_conditions.inst"

PRISMS_PF_END_NAMESPACE
