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
  set_dependencies_gradient_term_RHS(0, "grad(mu)");

  set_variable_name(1, "mu");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, AUXILIARY);
  set_dependencies_value_term_RHS(1, "c");
  set_dependencies_gradient_term_RHS(1, "grad(c)");

  set_variable_name(2, "f_tot");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);
  set_is_postprocessed_field(2, true);
  set_dependencies_value_term_RHS(2, "c, grad(c)");
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  scalarValue c   = variable_list.get_scalar_value(0);
  scalarGrad  mux = variable_list.get_scalar_gradient(1);

  scalarValue eq_c  = c;
  scalarGrad  eqx_c = -McV * this->user_inputs.temporal_discretization.dt * mux;

  variable_list.set_scalar_value_term(0, eq_c);
  variable_list.set_scalar_gradient_term(0, eqx_c);
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  if (this->current_index == 1)
    {
      scalarValue c  = variable_list.get_scalar_value(0);
      scalarGrad  cx = variable_list.get_scalar_gradient(0);

      scalarValue fcV = 4.0 * (c - 1.0) * (c - 0.5) * c;

      scalarValue eq_mu  = fcV;
      scalarGrad  eqx_mu = KcV * cx;

      variable_list.set_scalar_value_term(1, eq_mu);
      variable_list.set_scalar_gradient_term(1, eqx_mu);
    }
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_LHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_postprocess_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  scalarValue c  = variable_list.get_scalar_value(0);
  scalarGrad  cx = variable_list.get_scalar_gradient(0);

  scalarValue f_tot  = constV<number>(0.0);
  scalarValue f_chem = c * c * c * c - 2.0 * c * c * c + c * c;
  scalarValue f_grad = constV<number>(0.0);

  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += 0.5 * KcV * cx[i] * cx[j];
        }
    }
  f_tot = f_chem + f_grad;
  variable_list.set_scalar_value_term(2, f_tot);
}

INSTANTIATE_TRI_TEMPLATE(customPDE)

PRISMS_PF_END_NAMESPACE