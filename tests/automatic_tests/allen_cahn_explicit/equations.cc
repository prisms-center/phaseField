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
  set_variable_name(0, "n");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(0, "n");
  set_dependencies_gradient_term_RHS(0, "grad(n)");

  set_variable_name(1, "mg_n");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);
  set_is_postprocessed_field(1, true);
  set_dependencies_value_term_RHS(1, "grad(n)");

  set_variable_name(2, "f_tot");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);
  set_is_postprocessed_field(2, true);
  set_dependencies_value_term_RHS(2, "n, grad(n)");
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  scalarValue n  = variable_list.get_scalar_value(0);
  scalarGrad  nx = variable_list.get_scalar_gradient(0);

  scalarValue fnV   = 4.0 * n * (n - 1.0) * (n - 0.5);
  scalarValue eq_n  = n - user_inputs.temporal_discretization.dt * MnV * fnV;
  scalarGrad  eqx_n = -user_inputs.temporal_discretization.dt * KnV * MnV * nx;

  variable_list.set_scalar_value_term(0, eq_n);
  variable_list.set_scalar_gradient_term(0, eqx_n);
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

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
  scalarValue n  = variable_list.get_scalar_value(0);
  scalarGrad  nx = variable_list.get_scalar_gradient(0);

  scalarValue f_tot  = constV(static_cast<number>(0.0));
  scalarValue f_chem = n * n * n * n - 2.0 * n * n * n + n * n;
  scalarValue f_grad = constV(static_cast<number>(0.0));
  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += 0.5 * KnV * nx[i] * nx[j];
        }
    }
  f_tot = f_chem + f_grad;

  variable_list.set_scalar_value_term(1, std::sqrt(nx[0] * nx[0] + nx[1] * nx[1]));
  variable_list.set_scalar_value_term(2, f_tot);
}

INSTANTIATE_TRI_TEMPLATE(customPDE)

PRISMS_PF_END_NAMESPACE