// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
CustomAttributeLoader::load_variable_attributes()
{
  set_variable_name(0, "c");
  set_variable_type(0, Scalar);
  set_variable_equation_type(0, ImplicitTimeDependent);
  set_dependencies_value_term_rhs(0, "c, old_1(c)");
  set_dependencies_gradient_term_rhs(0, "grad(mu)");
  set_dependencies_value_term_lhs(0, "change(c)");
  set_dependencies_gradient_term_lhs(0, "");

  set_variable_name(1, "mu");
  set_variable_type(1, Scalar);
  set_variable_equation_type(1, Auxiliary);
  set_dependencies_value_term_rhs(1, "c");
  set_dependencies_gradient_term_rhs(1, "grad(c)");

  set_variable_name(2, "f_tot");
  set_variable_type(2, Scalar);
  set_variable_equation_type(2, ExplicitTimeDependent);
  set_is_postprocessed_field(2, true);
  set_dependencies_value_term_rhs(2, "c, grad(c)");
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
  [[maybe_unused]] Types::Index current_index) const
{
  if (current_index == 0)
    {
      ScalarValue c     = variable_list.get_scalar_value(0);
      ScalarValue old_c = variable_list.get_scalar_value(0, OldOne);
      ScalarGrad  mux   = variable_list.get_scalar_gradient(1);

      ScalarValue eq_c  = c - old_c;
      ScalarGrad  eqx_c = -McV * this->get_timestep() * mux;

      variable_list.set_scalar_value_term(0, eq_c);
      variable_list.set_scalar_gradient_term(0, eqx_c);
    }
  if (current_index == 1)
    {
      ScalarValue c  = variable_list.get_scalar_value(0);
      ScalarGrad  cx = variable_list.get_scalar_gradient(0);

      ScalarValue fcV = 4.0 * (c - 1.0) * (c - 0.5) * c;

      ScalarValue eq_mu  = fcV;
      ScalarGrad  eqx_mu = KcV * cx;

      variable_list.set_scalar_value_term(1, eq_mu);
      variable_list.set_scalar_gradient_term(1, eqx_mu);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index current_index) const
{
  if (current_index == 0)
    {
      ScalarValue change_c = variable_list.get_scalar_value(0, Change);

      variable_list.set_scalar_value_term(0, change_c, Change);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  ScalarValue c  = variable_list.get_scalar_value(0);
  ScalarGrad  cx = variable_list.get_scalar_gradient(0);

  ScalarValue f_tot  = 0.0;
  ScalarValue f_chem = c * c * c * c - 2.0 * c * c * c + c * c;
  ScalarValue f_grad = 0.0;

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_grad += 0.5 * KcV * cx[i] * cx[j];
        }
    }
  f_tot = f_chem + f_grad;
  variable_list.set_scalar_value_term(2, f_tot);
}

INSTANTIATE_TRI_TEMPLATE(CustomPDE)

PRISMS_PF_END_NAMESPACE
