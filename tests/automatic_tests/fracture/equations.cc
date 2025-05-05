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

      dealii::Tensor<2, CIJ_tensor_size, scalarValue> CIJ;
      for (unsigned int i = 0; i < CIJ_tensor_size; i++)
        {
          for (unsigned int j = 0; j < CIJ_tensor_size; j++)
            {
              CIJ[i][j] = CIJ_base[i][j] * Ex * (1.0 - 2.0 * n + n * n);
            }
        }
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

      dealii::Tensor<2, CIJ_tensor_size, scalarValue> CIJ;
      for (unsigned int i = 0; i < CIJ_tensor_size; i++)
        {
          for (unsigned int j = 0; j < CIJ_tensor_size; j++)
            {
              CIJ[i][j] = CIJ_base[i][j] * Ex;
            }
        }
      vectorGrad stress;
      compute_stress<dim, scalarValue>(CIJ, ux, stress);
      scalarValue elastic_energy = constV<number>(0.0);
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

      dealii::Tensor<2, CIJ_tensor_size, scalarValue> CIJ;
      for (unsigned int i = 0; i < CIJ_tensor_size; i++)
        {
          for (unsigned int j = 0; j < CIJ_tensor_size; j++)
            {
              CIJ[i][j] = CIJ_base[i][j] * Ex * (1.0 - 2.0 * n + n * n);
            }
        }
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
{}

INSTANTIATE_TRI_TEMPLATE(customPDE)

PRISMS_PF_END_NAMESPACE
