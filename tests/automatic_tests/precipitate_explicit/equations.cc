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
  set_variable_name(0, "c");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(0, "c");
  set_dependencies_gradient_term_RHS(
    0,
    "c, grad(c), n1, grad(n1), n2, grad(n2), n3, grad(n3), grad(u)");

  set_variable_name(1, "n1");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(1, "c, n1, n2, n3, grad(u)");
  set_dependencies_gradient_term_RHS(1, "grad(n1)");

  set_variable_name(2, "n2");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(2, "c, n1, n2, n3, grad(u)");
  set_dependencies_gradient_term_RHS(2, "grad(n2)");

  set_variable_name(3, "n3");
  set_variable_type(3, SCALAR);
  set_variable_equation_type(3, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(3, "c, n1, n2, n3, grad(u)");
  set_dependencies_gradient_term_RHS(3, "grad(n3)");

  set_variable_name(4, "u");
  set_variable_type(4, VECTOR);
  set_variable_equation_type(4, TIME_INDEPENDENT);

  set_dependencies_value_term_RHS(4, "");
  set_dependencies_gradient_term_RHS(4, "n1, n2, n3, grad(u)");
  set_dependencies_value_term_LHS(4, "");
  set_dependencies_gradient_term_LHS(4, "n1, n2, n3, grad(change(u))");

  set_variable_name(5, "f_tot");
  set_variable_type(5, SCALAR);
  set_variable_equation_type(5, EXPLICIT_TIME_DEPENDENT);

  set_dependencies_value_term_RHS(
    5,
    "c, grad(c), n1, grad(n1), n2, grad(n2), n3, grad(n3), grad(u)");
  set_dependencies_gradient_term_RHS(5, "");
  set_is_postprocessed_field(5, true);
}

#define compute_hV(n) (10.0 * n * n * n - 15.0 * n * n * n * n + 6.0 * n * n * n * n * n);
#define compute_hnV(n) (30.0 * n * n - 60.0 * n * n * n + 30.0 * n * n * n * n);

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  scalarValue c   = variable_list.get_scalar_value(0);
  scalarGrad  cx  = variable_list.get_scalar_gradient(0);
  scalarValue n1  = variable_list.get_scalar_value(1);
  scalarGrad  n1x = variable_list.get_scalar_gradient(1);
  scalarValue n2  = variable_list.get_scalar_value(2);
  scalarGrad  n2x = variable_list.get_scalar_gradient(2);
  scalarValue n3  = variable_list.get_scalar_value(3);
  scalarGrad  n3x = variable_list.get_scalar_gradient(3);
  vectorGrad  ux  = variable_list.get_vector_symmetric_gradient(4);

  // Free energy expressions and interpolation functions
  scalarValue faV   = A0 + A1 * c + A2 * c * c + A3 * c * c * c + A4 * c * c * c * c;
  scalarValue facV  = A1 + 2.0 * A2 * c + 3.0 * A3 * c * c + 4.0 * A4 * c * c * c;
  scalarValue faccV = 2.0 * A2 + 6.0 * A3 * c + 12.0 * A4 * c * c;
  scalarValue fbV   = B2 * c * c + B1 * c + B0;
  scalarValue fbcV  = 2.0 * B2 * c + B1;
  scalarValue fbccV = 2.0 * B2;
  scalarValue h1V   = compute_hV(n1);
  scalarValue h2V   = compute_hV(n2);
  scalarValue h3V   = compute_hV(n3);

  scalarValue hn1V = compute_hnV(n1);
  scalarValue hn2V = compute_hnV(n2);
  scalarValue hn3V = compute_hnV(n3);

  // Compute strain
  vectorGrad strain = ux - (sfts_const1 * h1V + sfts_const2 * h2V + sfts_const3 * h3V);

  // Compute stress
  vectorGrad stress;
  if (n_dependent_stiffness == true)
    {
      scalarValue                                     sum_hV = h1V + h2V + h3V;
      dealii::Tensor<2, CIJ_tensor_size, scalarValue> CIJ_combined =
        CIJ_Mg * (1.0 - sum_hV) + CIJ_Beta * sum_hV;
      compute_stress<dim, scalarValue>(CIJ_combined, strain, stress);
    }
  else
    {
      compute_stress<dim, scalarValue>(CIJ_Mg, strain, stress);
    }

  // Compute one of the stress terms in the order parameter chemical potential,
  // nDependentMisfitACp = C*(E-E0)*(E0_p*Hn)
  scalarValue nDependentMisfitAC1 = 0.0;
  scalarValue nDependentMisfitAC2 = 0.0;
  scalarValue nDependentMisfitAC3 = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          nDependentMisfitAC1 += stress[i][j] * sfts_const1[i][j];
          nDependentMisfitAC2 += stress[i][j] * sfts_const2[i][j];
          nDependentMisfitAC3 += stress[i][j] * sfts_const3[i][j];
        }
    }

  nDependentMisfitAC1 *= -hn1V;
  nDependentMisfitAC2 *= -hn2V;
  nDependentMisfitAC3 *= -hn3V;

  // Compute the other stress term in the order parameter chemical potential,
  // heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
  scalarValue heterMechAC1 = 0.0;
  scalarValue heterMechAC2 = 0.0;
  scalarValue heterMechAC3 = 0.0;
  if (n_dependent_stiffness == true)
    {
      vectorGrad stress_2;
      compute_stress<dim, scalarValue>(CIJ_Beta - CIJ_Mg, strain, stress_2);
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              heterMechAC1 += stress_2[i][j] * strain[i][j];
            }
        }
      // Aside from HnpV, heterMechAC1, heterMechAC2, and heterMechAC3 are equal
      heterMechAC2 = 0.5 * hn2V * heterMechAC1;
      heterMechAC3 = 0.5 * hn3V * heterMechAC1;
      heterMechAC1 = 0.5 * hn1V * heterMechAC1;
    }

  // compute K*nx
  scalarGrad Knx1, Knx2, Knx3;
  for (unsigned int a = 0; a < dim; a++)
    {
      Knx1[a] = 0.0;
      Knx2[a] = 0.0;
      Knx3[a] = 0.0;
      for (unsigned int b = 0; b < dim; b++)
        {
          Knx1[a] += Kn1[a][b] * n1x[b];
          Knx2[a] += Kn2[a][b] * n2x[b];
          Knx3[a] += Kn3[a][b] * n3x[b];
        }
    }

  scalarValue sum_hV = h1V + h2V + h3V;

  // The terms in the govering equations
  scalarValue eq_c       = c;
  scalarGrad  eqx_c_temp = cx * ((1.0 - sum_hV) * faccV + sum_hV * fbccV) +
                          n1x * ((fbcV - facV) * hn1V) + n2x * ((fbcV - facV) * hn2V) +
                          n3x * ((fbcV - facV) * hn3V);
  scalarGrad  eqx_c = -this->get_timestep() * McV * eqx_c_temp;
  scalarValue eq_n1 = n1 - this->get_timestep() * Mn1V *
                             ((fbV - faV) * hn1V + nDependentMisfitAC1 + heterMechAC1);
  scalarValue eq_n2 = n2 - this->get_timestep() * Mn2V *
                             ((fbV - faV) * hn2V + nDependentMisfitAC2 + heterMechAC2);
  scalarValue eq_n3 = n3 - this->get_timestep() * Mn3V *
                             ((fbV - faV) * hn3V + nDependentMisfitAC3 + heterMechAC3);
  scalarGrad eqx_n1 = -this->get_timestep() * Mn1V * Knx1;
  scalarGrad eqx_n2 = -this->get_timestep() * Mn2V * Knx2;
  scalarGrad eqx_n3 = -this->get_timestep() * Mn3V * Knx3;

  variable_list.set_scalar_value_term(0, eq_c);
  variable_list.set_scalar_gradient_term(0, eqx_c);
  variable_list.set_scalar_value_term(1, eq_n1);
  variable_list.set_scalar_gradient_term(1, eqx_n1);
  variable_list.set_scalar_value_term(2, eq_n2);
  variable_list.set_scalar_gradient_term(2, eqx_n2);
  variable_list.set_scalar_value_term(3, eq_n3);
  variable_list.set_scalar_gradient_term(3, eqx_n3);
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  if (current_index == 4)
    {
      scalarValue n1 = variable_list.get_scalar_value(1);
      scalarValue n2 = variable_list.get_scalar_value(2);
      scalarValue n3 = variable_list.get_scalar_value(3);
      vectorGrad  ux = variable_list.get_vector_symmetric_gradient(4);

      // Interpolation functions
      scalarValue h1V = compute_hV(n1);
      scalarValue h2V = compute_hV(n2);
      scalarValue h3V = compute_hV(n3);

      // Compute strain
      vectorGrad strain =
        ux - (sfts_const1 * h1V + sfts_const2 * h2V + sfts_const3 * h3V);

      // Compute stress
      vectorGrad stress;
      if (n_dependent_stiffness == true)
        {
          scalarValue                                     sum_hV = h1V + h2V + h3V;
          dealii::Tensor<2, CIJ_tensor_size, scalarValue> CIJ_combined =
            CIJ_Mg * (1.0 - sum_hV) + CIJ_Beta * sum_hV;
          compute_stress<dim, scalarValue>(CIJ_combined, strain, stress);
        }
      else
        {
          compute_stress<dim, scalarValue>(CIJ_Mg, strain, stress);
        }

      variable_list.set_vector_gradient_term(4, -stress);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_LHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  if (current_index == 4)
    {
      scalarValue n1        = variable_list.get_scalar_value(1);
      scalarValue n2        = variable_list.get_scalar_value(2);
      scalarValue n3        = variable_list.get_scalar_value(3);
      vectorGrad  change_ux = variable_list.get_vector_symmetric_gradient(4, CHANGE);

      // Interpolation functions
      scalarValue h1V = compute_hV(n1);
      scalarValue h2V = compute_hV(n2);
      scalarValue h3V = compute_hV(n3);

      // Compute strain
      vectorGrad strain = change_ux;

      // Compute stress
      vectorGrad stress;
      if (n_dependent_stiffness == true)
        {
          scalarValue                                     sum_hV = h1V + h2V + h3V;
          dealii::Tensor<2, CIJ_tensor_size, scalarValue> CIJ_combined =
            CIJ_Mg * (1.0 - sum_hV) + CIJ_Beta * sum_hV;
          compute_stress<dim, scalarValue>(CIJ_combined, strain, stress);
        }
      else
        {
          compute_stress<dim, scalarValue>(CIJ_Mg, strain, stress);
        }

      variable_list.set_vector_gradient_term(4, stress, CHANGE);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
customPDE<dim, degree, number>::compute_postprocess_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  scalarValue c   = variable_list.get_scalar_value(0);
  scalarValue n1  = variable_list.get_scalar_value(1);
  scalarGrad  n1x = variable_list.get_scalar_gradient(1);
  scalarValue n2  = variable_list.get_scalar_value(2);
  scalarGrad  n2x = variable_list.get_scalar_gradient(2);
  scalarValue n3  = variable_list.get_scalar_value(3);
  scalarGrad  n3x = variable_list.get_scalar_gradient(3);
  vectorGrad  ux  = variable_list.get_vector_symmetric_gradient(4);

  // Free energy expressions and interpolation functions
  scalarValue faV    = A0 + A1 * c + A2 * c * c + A3 * c * c * c + A4 * c * c * c * c;
  scalarValue fbV    = B2 * c * c + B1 * c + B0;
  scalarValue h1V    = compute_hV(n1);
  scalarValue h2V    = compute_hV(n2);
  scalarValue h3V    = compute_hV(n3);
  scalarValue sum_hV = h1V + h2V + h3V;

  scalarValue f_chem = (1.0 - sum_hV) * faV + sum_hV * fbV;

  scalarValue f_grad = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_grad += 0.5 * Kn1[i][j] * n1x[i] * n1x[j] +
                    0.5 * Kn2[i][j] * n2x[i] * n2x[j] + 0.5 * Kn3[i][j] * n3x[i] * n3x[j];
        }
    }

  // Compute strain
  vectorGrad strain = ux - (sfts_const1 * h1V + sfts_const2 * h2V + sfts_const3 * h3V);

  // Compute stress
  vectorGrad stress;
  if (n_dependent_stiffness == true)
    {
      dealii::Tensor<2, CIJ_tensor_size, scalarValue> CIJ_combined =
        CIJ_Mg * (1.0 - sum_hV) + CIJ_Beta * sum_hV;
      compute_stress<dim, scalarValue>(CIJ_combined, strain, stress);
    }
  else
    {
      compute_stress<dim, scalarValue>(CIJ_Mg, strain, stress);
    }

  scalarValue f_el = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_el += 0.5 * stress[i][j] * strain[i][j];
        }
    }

  scalarValue f_tot = f_chem + f_grad + f_el;

  variable_list.set_scalar_value_term(5, f_tot);
}

INSTANTIATE_TRI_TEMPLATE(customPDE)

PRISMS_PF_END_NAMESPACE
