// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
initialCondition<dim, degree>::initialCondition(
  const unsigned int                                            &_index,
  const fieldType                                               &field_type,
  const std::shared_ptr<const PDEOperator<dim, degree, double>> &_pde_operator)
  : dealii::Function<dim>((field_type == fieldType::VECTOR) ? dim : 1)
  , index(_index)
  , pde_operator(_pde_operator)
{}

// NOLINTBEGIN(readability-identifier-length)

template <unsigned int dim, unsigned int degree>
void
initialCondition<dim, degree>::vector_value(const dealii::Point<dim> &p,
                                            dealii::Vector<double>   &value) const
{
  // Initialize passed variables to zero
  dealii::Vector<double> vector_value(dim);

  // Pass variables to user-facing function to evaluate
  for (unsigned int i = 0; i < dim; i++)
    {
      pde_operator->set_initial_condition(index, i, p, vector_value[0], vector_value[i]);
    }

  value = vector_value;
}

// NOLINTEND(readability-identifier-length)

INSTANTIATE_BI_TEMPLATE(initialCondition)

PRISMS_PF_END_NAMESPACE