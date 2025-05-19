// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include "custom_pde.h"

#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

void
CustomAttributeLoader::loadVariableAttributes()
{
  set_variable_name(0, "n");
  set_variable_type(0, Scalar);
  set_variable_equation_type(0, ImplicitTimeDependent);
  set_dependencies_value_term_rhs(0, "n, old_1(n)");
  set_dependencies_gradient_term_rhs(0, "grad(n)");
  set_dependencies_value_term_lhs(0, "change(n)");
  set_dependencies_gradient_term_lhs(0, "grad(change(n))");

  set_variable_name(1, "mg_n");
  set_variable_type(1, Scalar);
  set_variable_equation_type(1, ExplicitTimeDependent);
  set_is_postprocessed_field(1, true);
  set_dependencies_value_term_rhs(1, "grad(n)");

  set_variable_name(2, "f_tot");
  set_variable_type(2, Scalar);
  set_variable_equation_type(2, ExplicitTimeDependent);
  set_is_postprocessed_field(2, true);
  set_dependencies_value_term_rhs(2, "n, grad(n)");
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  if (current_index == 0)
    {
      scalarValue n     = variable_list.get_scalar_value(0);
      scalarValue old_n = variable_list.get_scalar_value(0, OldOne);
      scalarGrad  nx    = variable_list.get_scalar_gradient(0);

      scalarValue fnV   = 4.0 * n * (n - 1.0) * (n - 0.5);
      scalarValue eq_n  = old_n - n - this->get_timestep() * MnV * fnV;
      scalarGrad  eqx_n = -this->get_timestep() * KnV * MnV * nx;

      variable_list.set_scalar_value_term(0, eq_n);
      variable_list.set_scalar_gradient_term(0, eqx_n);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  if (current_index == 0)
    {
      scalarValue change_n  = variable_list.get_scalar_value(0, Change);
      scalarGrad  change_nx = variable_list.get_scalar_gradient(0, Change);

      scalarValue fnV          = 4.0 * change_n * (change_n - 1.0) * (change_n - 0.5);
      scalarValue eq_change_n  = change_n + this->get_timestep() * MnV * fnV;
      scalarGrad  eqx_change_n = this->get_timestep() * KnV * MnV * change_nx;

      variable_list.set_scalar_value_term(0, eq_change_n, Change);
      variable_list.set_scalar_gradient_term(0, eqx_change_n, Change);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  scalarValue n  = variable_list.get_scalar_value(0);
  scalarGrad  nx = variable_list.get_scalar_gradient(0);

  scalarValue f_tot  = 0.0;
  scalarValue f_chem = n * n * n * n - 2.0 * n * n * n + n * n;
  scalarValue f_grad = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
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
