// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/config.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>

PRISMS_PF_BEGIN_NAMESPACE

void
customAttributeLoader::loadVariableAttributes()
{
  set_variable_name(0, "c");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(0, "c");
  set_dependencies_gradient_term_RHS(0, "");

  set_variable_name(1, "n1");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(1, "n1");
  set_dependencies_gradient_term_RHS(1, "");

  set_variable_name(2, "n2");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(2, "n2");
  set_dependencies_gradient_term_RHS(2, "");

  set_variable_name(3, "n3");
  set_variable_type(3, SCALAR);
  set_variable_equation_type(3, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(3, "n3");
  set_dependencies_gradient_term_RHS(3, "");

  set_variable_name(4, "u");
  set_variable_type(4, VECTOR);
  set_variable_equation_type(4, TIME_INDEPENDENT);
  set_dependencies_value_term_RHS(4, "");
  set_dependencies_gradient_term_RHS(4, "grad(u)");
  set_dependencies_value_term_LHS(4, "");
  set_dependencies_gradient_term_LHS(4, "grad(change(u))");
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  scalarValue c  = variable_list.get_scalar_value(0);
  scalarValue n1 = variable_list.get_scalar_value(1);
  scalarValue n2 = variable_list.get_scalar_value(2);
  scalarValue n3 = variable_list.get_scalar_value(3);

  variable_list.set_scalar_value_term(0, c);
  variable_list.set_scalar_value_term(1, n1);
  variable_list.set_scalar_value_term(2, n2);
  variable_list.set_scalar_value_term(3, n3);
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  if (this->current_index == 4)
    {
      vectorGrad grad_u = variable_list.get_vector_gradient(4);

      variable_list.set_vector_gradient_term(4, -grad_u);
    }
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_LHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  if (this->current_index == 4)
    {
      vectorGrad change_grad_u = variable_list.get_vector_gradient(4, CHANGE);

      variable_list.set_vector_gradient_term(4, change_grad_u, CHANGE);
    }
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_postprocess_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

INSTANTIATE_TRI_TEMPLATE(customPDE)

PRISMS_PF_END_NAMESPACE
