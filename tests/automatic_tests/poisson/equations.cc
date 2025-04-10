// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
customAttributeLoader::loadVariableAttributes()
{
  set_variable_name(0, "u");
  set_variable_type(0, VECTOR);
  set_variable_equation_type(0, TIME_INDEPENDENT);
  set_dependencies_gradient_term_RHS(0, "grad(u)");
  set_dependencies_gradient_term_LHS(0, "grad(change(u))");
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  if (current_index == 0)
    {
      vectorGrad ux = variable_list.get_vector_gradient(0);

      variable_list.set_vector_gradient_term(0, -ux);
    }
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_LHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  if (current_index == 0)
    {
      vectorGrad change_ux = variable_list.get_vector_gradient(0, CHANGE);

      variable_list.set_vector_gradient_term(0, change_ux, CHANGE);
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
