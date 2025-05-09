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
  set_variable_name(0, "n");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "n, dndt");
  set_dependencies_gradient_term_RHS(0, "");

  set_variable_name(1, "u");
  set_variable_type(1, VECTOR);
  set_variable_equation_type(1, TIME_INDEPENDENT);

  set_dependencies_value_term_RHS(1, "");
  set_dependencies_gradient_term_RHS(1, "n, grad(u), Ex");
  set_dependencies_value_term_LHS(1, "");
  set_dependencies_gradient_term_LHS(1, "n, grad(change(u)), Ex");

  set_variable_name(2, "dndt");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, AUXILIARY);

  set_dependencies_value_term_RHS(2, "n, grad(u), Ex, Gx");
  set_dependencies_gradient_term_RHS(2, "grad(n), Gx");

  set_variable_name(3, "Ex");
  set_variable_type(3, SCALAR);
  set_variable_equation_type(3, CONSTANT);

  set_variable_name(4, "Gx");
  set_variable_type(4, SCALAR);
  set_variable_equation_type(4, CONSTANT);

  set_variable_name(5, "f_tot");
  set_variable_type(5, SCALAR);
  set_variable_equation_type(5, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(5, "n, grad(n), grad(u), Ex, Gx");
  set_dependencies_gradient_term_RHS(5, "");
  set_is_postprocessed_field(5, true);

  set_variable_name(6, "s11");
  set_variable_type(6, SCALAR);
  set_variable_equation_type(6, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(6, "n, grad(u), Ex");
  set_dependencies_gradient_term_RHS(6, "");
  set_is_postprocessed_field(6, true);

  set_variable_name(7, "s12");
  set_variable_type(7, SCALAR);
  set_variable_equation_type(7, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(7, "n, grad(u), Ex");
  set_dependencies_gradient_term_RHS(7, "");
  set_is_postprocessed_field(7, true);

  set_variable_name(8, "s22");
  set_variable_type(8, SCALAR);
  set_variable_equation_type(8, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(8, "n, grad(u), Ex");
  set_dependencies_gradient_term_RHS(8, "");
  set_is_postprocessed_field(8, true);

  set_variable_name(9, "e22");
  set_variable_type(9, SCALAR);
  set_variable_equation_type(9, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(9, "grad(u), Ex");
  set_dependencies_gradient_term_RHS(9, "");
  set_is_postprocessed_field(9, true);

  set_variable_name(10, "f_int");
  set_variable_type(10, SCALAR);
  set_variable_equation_type(10, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(10, "n, grad(n), grad(u), Ex, Gx");
  set_dependencies_gradient_term_RHS(10, "");
  set_is_postprocessed_field(10, true);

  set_variable_name(11, "f_el");
  set_variable_type(11, SCALAR);
  set_variable_equation_type(11, EXPLICIT_TIME_DEPENDENT);
  set_dependencies_value_term_RHS(11, "n, grad(u), Ex");
  set_dependencies_gradient_term_RHS(11, "");
  set_is_postprocessed_field(11, true);
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  scalarValue n    = variable_list.get_scalar_value(0);
  scalarValue dndt = variable_list.get_scalar_value(2);

  for (unsigned int j = 0; j < dndt.size(); ++j)
    {
      if (dndt[j] > 0.0)
        {
          dndt[j] = 0.0;
        }
      if (n[j] - dndt[j] * this->get_timestep() > 1.0)
        {
          dndt[j] = (n[j] - 1.0) / this->get_timestep();
        }
    }
  scalarValue eq_n = n - this->get_timestep() * Mn * dndt;

  variable_list.set_scalar_value_term(0, eq_n);
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  if (current_index == 1)
    {
      scalarValue n  = variable_list.get_scalar_value(0);
      vectorGrad  ux = variable_list.get_vector_symmetric_gradient(1);
      scalarValue Ex = variable_list.get_scalar_value(3);

      dealii::Tensor<2, voigt_tensor_size<dim>, scalarValue> CIJ =
        CIJ_base * Ex * (1.0 - 2.0 * n + n * n);
      vectorGrad stress;
      compute_stress<dim, scalarValue>(CIJ, ux, stress);

      variable_list.set_vector_gradient_term(1, -stress);
    }
  if (current_index == 2)
    {
      scalarValue n  = variable_list.get_scalar_value(0);
      scalarGrad  nx = variable_list.get_scalar_gradient(0);
      vectorGrad  ux = variable_list.get_vector_symmetric_gradient(1);
      scalarValue Ex = variable_list.get_scalar_value(3);
      scalarValue Gx = variable_list.get_scalar_value(4);

      dealii::Tensor<2, voigt_tensor_size<dim>, scalarValue> CIJ = CIJ_base * Ex;
      vectorGrad                                             stress;
      compute_stress<dim, scalarValue>(CIJ, ux, stress);
      scalarValue elastic_energy = 0.0;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              elastic_energy += 0.5 * stress[i][j] * ux[i][j];
            }
        }
      scalarValue eq_n =
        (2.0 * (n - 1.0) * elastic_energy + Gc0 * Gx * 3.0 / 8.0 / ell) * Mn;
      scalarGrad eqx_n = ell * nx * Gc0 * Gx * 3.0 / 8.0 * Mn;

      variable_list.set_scalar_value_term(2, eq_n);
      variable_list.set_scalar_gradient_term(2, eqx_n);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_LHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  if (current_index == 1)
    {
      scalarValue n         = variable_list.get_scalar_value(0);
      vectorGrad  ux_change = variable_list.get_vector_symmetric_gradient(1, CHANGE);
      scalarValue Ex        = variable_list.get_scalar_value(3);

      dealii::Tensor<2, voigt_tensor_size<dim>, scalarValue> CIJ =
        CIJ_base * Ex * (1.0 - 2.0 * n + n * n);
      vectorGrad stress;
      compute_stress<dim, scalarValue>(CIJ, ux_change, stress);

      variable_list.set_vector_gradient_term(1, stress, CHANGE);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_postprocess_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  scalarValue n  = variable_list.get_scalar_value(0);
  scalarGrad  nx = variable_list.get_scalar_gradient(0);
  vectorGrad  ux = variable_list.get_vector_symmetric_gradient(1);
  scalarValue Ex = variable_list.get_scalar_value(3);
  scalarValue Gx = variable_list.get_scalar_value(4);

  scalarValue f_int = Gc0 * n * Gx * 3.0 / 8.0 / ell + Gc0 * Gx * 0.5 * ell * nx * nx;

  dealii::Tensor<2, voigt_tensor_size<dim>, scalarValue> CIJ =
    CIJ_base * Ex * (1.0 - 2.0 * n + n * n);
  vectorGrad stress;
  compute_stress<dim, scalarValue>(CIJ, ux, stress);

  scalarValue f_el = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_el += 0.5 * stress[i][j] * ux[i][j];
        }
    }

  scalarValue f_tot = f_el + f_int;
  scalarValue s11   = stress[0][0];
  scalarValue s12   = stress[0][1];
  scalarValue s22   = stress[1][1];
  scalarValue e22   = ux[1][1];

  variable_list.set_scalar_value_term(5, f_tot);
  variable_list.set_scalar_value_term(6, s11);
  variable_list.set_scalar_value_term(7, s12);
  variable_list.set_scalar_value_term(8, s22);
  variable_list.set_scalar_value_term(9, e22);
  variable_list.set_scalar_value_term(10, f_int);
  variable_list.set_scalar_value_term(11, f_el);
}

INSTANTIATE_TRI_TEMPLATE(customPDE)

PRISMS_PF_END_NAMESPACE
