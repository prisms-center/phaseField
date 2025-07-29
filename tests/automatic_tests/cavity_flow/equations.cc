// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
CustomAttributeLoader::load_variable_attributes()
{
  set_variable_name(0, "u");
  set_variable_type(0, Vector);
  set_variable_equation_type(0, ExplicitTimeDependent);

  set_dependencies_value_term_rhs(0, "u, grad(P)");
  set_dependencies_gradient_term_rhs(0, "");

  set_variable_name(1, "u_star");
  set_variable_type(1, Vector);
  set_variable_equation_type(1, ExplicitTimeDependent);

  set_dependencies_value_term_rhs(1, "u, grad(u)");
  set_dependencies_gradient_term_rhs(1, "grad(u)");

  set_variable_name(2, "P");
  set_variable_type(2, Scalar);
  set_variable_equation_type(2, TimeIndependent);

  set_dependencies_value_term_rhs(2, "div(u_star)");
  set_dependencies_gradient_term_rhs(2, "grad(P)");
  set_dependencies_value_term_lhs(2, "");
  set_dependencies_gradient_term_lhs(2, "grad(change(P))");
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index                                               index) const
{
  if (index == 2)
    {
      ScalarValue div_u_star = variable_list.template get_divergence<ScalarValue>(1);
      ScalarGrad  Px         = variable_list.template get_gradient<ScalarGrad>(2);

      variable_list.set_value_term(2, div_u_star);
      variable_list.set_gradient_term(2, -Px);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index                                               index) const
{
  if (index == 2)
    {
      ScalarGrad change_Px = variable_list.template get_gradient<ScalarGrad>(2, Change);

      variable_list.set_gradient_term(2, change_Px);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

INSTANTIATE_TRI_TEMPLATE(CustomPDE)

PRISMS_PF_END_NAMESPACE
