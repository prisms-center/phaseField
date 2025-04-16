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
    "c, grad(c), n1, grad(n1), n2, grad(n2), n3, grad(n3), grad(u), hess(u)");

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
  set_dependencies_gradient_term_RHS(4, "c, n1, n2, n3, grad(u)");
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

template <int dim, int degree, typename number>
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
  vectorGrad  ux  = variable_list.get_vector_gradient(4);

  vectorHess uxx;
  bool       c_dependent_misfit = false;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          if (std::abs(sfts_linear1[i][j]) > 1.0e-12)
            {
              c_dependent_misfit = true;
            }
        }
    }
  if (c_dependent_misfit == true)
    {
      uxx = variable_list.get_vector_hessian(4);
    }

  // Free energy expressions and interpolation functions
  scalarValue faV   = A0 + A1 * c + A2 * c * c + A3 * c * c * c + A4 * c * c * c * c;
  scalarValue facV  = A1 + 2.0 * A2 * c + 3.0 * A3 * c * c + 4.0 * A4 * c * c * c;
  scalarValue faccV = 2.0 * A2 + 6.0 * A3 * c + 12.0 * A4 * c * c;
  scalarValue fbV   = B2 * c * c + B1 * c + B0;
  scalarValue fbcV  = 2.0 * B2 * c + B1;
  scalarValue fbccV = 2.0 * B2;
  scalarValue h1V =
    10.0 * n1 * n1 * n1 - 15.0 * n1 * n1 * n1 * n1 + 6.0 * n1 * n1 * n1 * n1 * n1;
  scalarValue h2V =
    10.0 * n2 * n2 * n2 - 15.0 * n2 * n2 * n2 * n2 + 6.0 * n2 * n2 * n2 * n2 * n2;
  scalarValue h3V =
    10.0 * n3 * n3 * n3 - 15.0 * n3 * n3 * n3 * n3 + 6.0 * n3 * n3 * n3 * n3 * n3;
  scalarValue hn1V = 30.0 * n1 * n1 - 60.0 * n1 * n1 * n1 + 30.0 * n1 * n1 * n1 * n1;
  scalarValue hn2V = 30.0 * n2 * n2 - 60.0 * n2 * n2 * n2 + 30.0 * n2 * n2 * n2 * n2;
  scalarValue hn3V = 30.0 * n3 * n3 - 60.0 * n3 * n3 * n3 + 30.0 * n3 * n3 * n3 * n3;

  // Calculate the stress-free transformation strain and its derivatives at the
  // quadrature point
  vectorGrad sfts1, sfts1c, sfts1cc, sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c + b_p
          sfts1[i][j]   = sfts_linear1[i][j] * c + sfts_const1[i][j];
          sfts1c[i][j]  = sfts_linear1[i][j];
          sfts1cc[i][j] = constV<number>(0.0);

          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c + b_p
          sfts2[i][j]   = sfts_linear2[i][j] * c + sfts_const2[i][j];
          sfts2c[i][j]  = sfts_linear1[i][j];
          sfts2cc[i][j] = constV<number>(0.0);

          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c + b_p
          sfts3[i][j]   = sfts_linear3[i][j] * c + sfts_const3[i][j];
          sfts3c[i][j]  = sfts_linear3[i][j];
          sfts3cc[i][j] = constV<number>(0.0);
        }
    }

  // compute strain_2=(E-E0)
  scalarValue strain_2[dim][dim], stress[dim][dim];

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          strain_2[i][j] = 0.5 * (ux[i][j] + ux[j][i]) -
                           (sfts1[i][j] * h1V + sfts2[i][j] * h2V + sfts3[i][j] * h3V);
        }
    }

  // compute stress
  // stress=C*(E-E0)
  //  Compute stress tensor (which is equal to the residual, Rux)
  scalarValue CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

  if (n_dependent_stiffness == true)
    {
      scalarValue sum_hV;
      sum_hV = h1V + h2V + h3V;
      for (unsigned int i = 0; i < 2 * dim - 1 + dim / 3; i++)
        {
          for (unsigned int j = 0; j < 2 * dim - 1 + dim / 3; j++)
            {
              CIJ_combined[i][j] =
                CIJ_Mg[i][j] * (1.0 - sum_hV) + CIJ_Beta[i][j] * sum_hV;
            }
        }
      compute_stress<dim, scalarValue>(CIJ_combined, strain_2, stress);
    }
  else
    {
      compute_stress<dim, scalarValue>(CIJ_Mg, strain_2, stress);
    }

  // Compute one of the stress terms in the order parameter chemical potential,
  // nDependentMisfitACp = C*(E-E0)*(E0_p*Hn)
  scalarValue nDependentMisfitAC1 = constV<number>(0.0);
  scalarValue nDependentMisfitAC2 = constV<number>(0.0);
  scalarValue nDependentMisfitAC3 = constV<number>(0.0);

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          nDependentMisfitAC1 += stress[i][j] * sfts1[i][j];
          nDependentMisfitAC2 += stress[i][j] * sfts2[i][j];
          nDependentMisfitAC3 += stress[i][j] * sfts3[i][j];
        }
    }

  nDependentMisfitAC1 *= -hn1V;
  nDependentMisfitAC2 *= -hn2V;
  nDependentMisfitAC3 *= -hn3V;

  // Compute the other stress term in the order parameter chemical potential,
  // heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
  scalarValue heterMechAC1 = constV<number>(0.0);
  scalarValue heterMechAC2 = constV<number>(0.0);
  scalarValue heterMechAC3 = constV<number>(0.0);
  scalarValue stress_2[dim][dim];

  if (n_dependent_stiffness == true)
    {
      // compute_stress<dim,number>(CIJ_diff, strain_2, stress_2);
      compute_stress<dim, scalarValue>(CIJ_Beta - CIJ_Mg, strain_2, stress_2);
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              heterMechAC1 += stress_2[i][j] * strain_2[i][j];
            }
        }
      // Aside from HnpV, heterMechAC1, heterMechAC2, and heterMechAC3 are equal
      heterMechAC2 = 0.5 * hn2V * heterMechAC1;
      heterMechAC3 = 0.5 * hn3V * heterMechAC1;

      heterMechAC1 = 0.5 * hn1V * heterMechAC1;
    }

  // compute the stress term in the gradient of the concentration chemical
  // potential, grad_mu_el = [C*(E-E0)*E0c]x, must be a vector with length dim
  scalarGrad grad_mu_el;

  if (c_dependent_misfit == true)
    {
      scalarValue strain_3[dim][dim], stress_3[dim][dim];

      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              strain_3[i][j] =
                -(sfts1c[i][j] * h1V + sfts2c[i][j] * h2V + sfts3c[i][j] * h3V);
            }
        }

      if (n_dependent_stiffness == true)
        {
          compute_stress<dim, scalarValue>(CIJ_combined, strain_3, stress_3);
        }
      else
        {
          compute_stress<dim, scalarValue>(CIJ_Mg, strain_3, stress_3);
        }

      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              for (unsigned int k = 0; k < dim; k++)
                {
                  grad_mu_el[k] +=
                    stress_3[i][j] *
                    (0.5 * (uxx[i][j][k] + uxx[j][i][k]) + strain_3[i][j] * cx[k] -
                     (sfts1[i][j] * hn1V * n1x[k] + sfts2[i][j] * hn2V * n2x[k] +
                      sfts3[i][j] * hn3V * n3x[k]));

                  grad_mu_el[k] +=
                    -stress[i][j] *
                    (sfts1c[i][j] * hn1V * n1x[k] + sfts2c[i][j] * hn2V * n2x[k] +
                     sfts3c[i][j] * hn3V * n3x[k] +
                     (sfts1cc[i][j] * h1V + sfts2cc[i][j] * h2V + sfts3cc[i][j] * h3V) *
                       cx[k]);

                  if (n_dependent_stiffness == true)
                    {
                      grad_mu_el[k] += stress_2[i][j] * strain_3[i][j] *
                                       (hn1V * n1x[k] + hn2V * n2x[k] + hn3V * n3x[k]);
                    }
                }
            }
        }
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

  // The terms in the govering equations
  scalarValue eq_c = (c);
  scalarGrad  eqx_c_temp =
    (cx * ((1.0 - h1V - h2V - h3V) * faccV + (h1V + h2V + h3V) * fbccV) +
     n1x * ((fbcV - facV) * hn1V) + n2x * ((fbcV - facV) * hn2V) +
     n3x * ((fbcV - facV) * hn3V) + grad_mu_el);
  scalarGrad  eqx_c  = (-this->get_timestep() * McV * eqx_c_temp);
  scalarValue eq_n1  = (n1 - this->get_timestep() * Mn1V *
                              ((fbV - faV) * hn1V + nDependentMisfitAC1 + heterMechAC1));
  scalarValue eq_n2  = (n2 - this->get_timestep() * Mn2V *
                              ((fbV - faV) * hn2V + nDependentMisfitAC2 + heterMechAC2));
  scalarValue eq_n3  = (n3 - this->get_timestep() * Mn3V *
                              ((fbV - faV) * hn3V + nDependentMisfitAC3 + heterMechAC3));
  scalarGrad  eqx_n1 = (-this->get_timestep() * Mn1V * Knx1);
  scalarGrad  eqx_n2 = (-this->get_timestep() * Mn2V * Knx2);
  scalarGrad  eqx_n3 = (-this->get_timestep() * Mn3V * Knx3);

  variable_list.set_scalar_value_term(0, eq_c);
  variable_list.set_scalar_gradient_term(0, eqx_c);
  variable_list.set_scalar_value_term(1, eq_n1);
  variable_list.set_scalar_gradient_term(1, eqx_n1);
  variable_list.set_scalar_value_term(2, eq_n2);
  variable_list.set_scalar_gradient_term(2, eqx_n2);
  variable_list.set_scalar_value_term(3, eq_n3);
  variable_list.set_scalar_gradient_term(3, eqx_n3);
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  if (current_index == 4)
    {
      scalarValue c  = variable_list.get_scalar_value(0);
      scalarValue n1 = variable_list.get_scalar_value(1);
      scalarValue n2 = variable_list.get_scalar_value(2);
      scalarValue n3 = variable_list.get_scalar_value(3);
      vectorGrad  ux = variable_list.get_vector_gradient(4);

      // Interpolation functions
      scalarValue h1V =
        10.0 * n1 * n1 * n1 - 15.0 * n1 * n1 * n1 * n1 + 6.0 * n1 * n1 * n1 * n1 * n1;
      scalarValue h2V =
        10.0 * n2 * n2 * n2 - 15.0 * n2 * n2 * n2 * n2 + 6.0 * n2 * n2 * n2 * n2 * n2;
      scalarValue h3V =
        10.0 * n3 * n3 * n3 - 15.0 * n3 * n3 * n3 * n3 + 6.0 * n3 * n3 * n3 * n3 * n3;

      // Calculate the stress-free transformation strain and its derivatives at the
      // quadrature point
      vectorGrad sfts1, sfts2, sfts3;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              // Polynomial fits for the stress-free transformation strains, of the
              // form: sfts = a_p * c + b_p
              sfts1[i][j] = sfts_linear1[i][j] * c + sfts_const1[i][j];
              // Polynomial fits for the stress-free transformation strains, of the
              // form: sfts = a_p * c + b_p
              sfts2[i][j] = sfts_linear2[i][j] * c + sfts_const2[i][j];
              // Polynomial fits for the stress-free transformation strains, of the
              // form: sfts = a_p * c + b_p
              sfts3[i][j] = sfts_linear3[i][j] * c + sfts_const3[i][j];
            }
        }

      // compute strain_2=(E-E0)
      scalarValue strain_2[dim][dim], stress[dim][dim];

      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              strain_2[i][j] =
                0.5 * (ux[i][j] + ux[j][i]) -
                (sfts1[i][j] * h1V + sfts2[i][j] * h2V + sfts3[i][j] * h3V);
            }
        }

      // compute stress
      // stress=C*(E-E0)
      //  Compute stress tensor (which is equal to the residual, Rux)
      scalarValue CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

      if (n_dependent_stiffness == true)
        {
          scalarValue sum_hV;
          sum_hV = h1V + h2V + h3V;
          for (unsigned int i = 0; i < CIJ_tensor_size; i++)
            {
              for (unsigned int j = 0; j < CIJ_tensor_size; j++)
                {
                  CIJ_combined[i][j] =
                    CIJ_Mg[i][j] * (1.0 - sum_hV) + CIJ_Beta[i][j] * sum_hV;
                }
            }
          compute_stress<dim, scalarValue>(CIJ_combined, strain_2, stress);
        }
      else
        {
          compute_stress<dim, scalarValue>(CIJ_Mg, strain_2, stress);
        }

      vectorGrad eqx_u;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              eqx_u[i][j] = -stress[i][j];
            }
        }

      variable_list.set_vector_gradient_term(4, eqx_u);
    }
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_LHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] types::index current_index) const
{
  if (current_index == 4)
    {
      scalarValue n1  = variable_list.get_scalar_value(1);
      scalarValue n2  = variable_list.get_scalar_value(2);
      scalarValue n3  = variable_list.get_scalar_value(3);
      vectorGrad  Dux = variable_list.get_vector_gradient(4, CHANGE);

      vectorGrad eqx_Du;

      // Interpolation functions

      scalarValue h1V =
        10.0 * n1 * n1 * n1 - 15.0 * n1 * n1 * n1 * n1 + 6.0 * n1 * n1 * n1 * n1 * n1;
      scalarValue h2V =
        10.0 * n2 * n2 * n2 - 15.0 * n2 * n2 * n2 * n2 + 6.0 * n2 * n2 * n2 * n2 * n2;
      scalarValue h3V =
        10.0 * n3 * n3 * n3 - 15.0 * n3 * n3 * n3 * n3 + 6.0 * n3 * n3 * n3 * n3 * n3;

      // Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the
      // dealii "symmetrize" function
      vectorGrad E;
      E = dealii::symmetrize(Dux);

      // Compute stress tensor (which is equal to the residual, Rux)
      if (n_dependent_stiffness == true)
        {
          dealii::Tensor<2, CIJ_tensor_size, scalarValue> CIJ_combined;
          CIJ_combined = CIJ_Mg * (1.0 - h1V - h2V - h3V) + CIJ_Beta * (h1V + h2V + h3V);

          compute_stress<dim, scalarValue>(CIJ_combined, E, eqx_Du);
        }
      else
        {
          compute_stress<dim, scalarValue>(CIJ_Mg, E, eqx_Du);
        }

      variable_list.set_vector_gradient_term(4, eqx_Du, CHANGE);
    }
}

template <int dim, int degree, typename number>
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
  vectorGrad  ux  = variable_list.get_vector_gradient(4);

  scalarValue f_tot = constV<number>(0.0);

  // Free energy expressions and interpolation functions
  scalarValue faV = A0 + A1 * c + A2 * c * c + A3 * c * c * c + A4 * c * c * c * c;
  scalarValue fbV = B2 * c * c + B1 * c + B0;
  scalarValue h1V =
    10.0 * n1 * n1 * n1 - 15.0 * n1 * n1 * n1 * n1 + 6.0 * n1 * n1 * n1 * n1 * n1;
  scalarValue h2V =
    10.0 * n2 * n2 * n2 - 15.0 * n2 * n2 * n2 * n2 + 6.0 * n2 * n2 * n2 * n2 * n2;
  scalarValue h3V =
    10.0 * n3 * n3 * n3 - 15.0 * n3 * n3 * n3 * n3 + 6.0 * n3 * n3 * n3 * n3 * n3;

  scalarValue f_chem = (1.0 - (h1V + h2V + h3V)) * faV + (h1V + h2V + h3V) * fbV;

  scalarValue f_grad = constV<number>(0.0);

  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += 0.5 * Kn1[i][j] * n1x[i] * n1x[j];
        }
    }

  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += 0.5 * Kn2[i][j] * n2x[i] * n2x[j];
        }
    }

  for (int i = 0; i < dim; i++)
    {
      for (int j = 0; j < dim; j++)
        {
          f_grad += 0.5 * Kn3[i][j] * n3x[i] * n3x[j];
        }
    }

  // Calculate the stress-free transformation strain and its derivatives at the
  // quadrature point
  vectorGrad sfts1, sfts1c, sfts1cc, sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c + b_p
          sfts1[i][j]   = sfts_linear1[i][j] * c + sfts_const1[i][j];
          sfts1c[i][j]  = sfts_linear1[i][j];
          sfts1cc[i][j] = constV<number>(0.0);

          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c + b_p
          sfts2[i][j]   = sfts_linear2[i][j] * c + sfts_const2[i][j];
          sfts2c[i][j]  = sfts_linear1[i][j];
          sfts2cc[i][j] = constV<number>(0.0);

          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c + b_p
          sfts3[i][j]   = sfts_linear3[i][j] * c + sfts_const3[i][j];
          sfts3c[i][j]  = sfts_linear3[i][j];
          sfts3cc[i][j] = constV<number>(0.0);
        }
    }

  // compute strain_2=(E-E0)
  scalarValue strain_2[dim][dim], stress[dim][dim];

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          // strain_2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V +
          // sfts2[i][j]*h2V + sfts3[i][j]*h3V);
          strain_2[i][j] = 0.5 * (ux[i][j] + ux[j][i]) -
                           (sfts1[i][j] * h1V + sfts2[i][j] * h2V + sfts3[i][j] * h3V);
        }
    }

  // compute stress
  // stress=C*(E-E0)
  scalarValue CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

  if (n_dependent_stiffness == true)
    {
      scalarValue sum_hV;
      sum_hV = h1V + h2V + h3V;
      for (unsigned int i = 0; i < 2 * dim - 1 + dim / 3; i++)
        {
          for (unsigned int j = 0; j < 2 * dim - 1 + dim / 3; j++)
            {
              CIJ_combined[i][j] =
                CIJ_Mg[i][j] * (1.0 - sum_hV) + CIJ_Beta[i][j] * sum_hV;
            }
        }
      compute_stress<dim, scalarValue>(CIJ_combined, strain_2, stress);
    }
  else
    {
      compute_stress<dim, scalarValue>(CIJ_Mg, strain_2, stress);
    }

  scalarValue f_el = constV<number>(0.0);

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_el += 0.5 * stress[i][j] * strain_2[i][j];
        }
    }

  f_tot = f_chem + f_grad + f_el;

  variable_list.set_scalar_value_term(5, f_tot);
}

INSTANTIATE_TRI_TEMPLATE(customPDE)

PRISMS_PF_END_NAMESPACE
