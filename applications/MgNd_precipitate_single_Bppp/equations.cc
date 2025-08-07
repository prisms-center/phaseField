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
  set_variable_name(0, "c");
  set_variable_type(0, Scalar);
  set_variable_equation_type(0, ExplicitTimeDependent);

  set_dependencies_value_term_rhs(0, "c");
  set_dependencies_gradient_term_rhs(0, "n1, grad(mu)");

  set_variable_name(1, "mu");
  set_variable_type(1, Scalar);
  set_variable_equation_type(1, Auxiliary);

  set_dependencies_value_term_rhs(1, "c, n1, grad(u)");
  set_dependencies_gradient_term_rhs(1, "");

  set_variable_name(2, "n1");
  set_variable_type(2, Scalar);
  set_variable_equation_type(2, ExplicitTimeDependent);

  set_dependencies_value_term_rhs(2, "c, n1, grad(u)");
  set_dependencies_gradient_term_rhs(2, "grad(n1)");

  set_variable_name(3, "u");
  set_variable_type(3, Vector);
  set_variable_equation_type(3, ExplicitTimeDependent);

  set_dependencies_value_term_rhs(3, "");
  set_dependencies_gradient_term_rhs(3, "n1, grad(u)");
  set_dependencies_value_term_lhs(3, "");
  set_dependencies_gradient_term_lhs(3, "n1, grad(change(u))");
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  ScalarValue c       = variable_list.template get_value<ScalarValue>(0);
  ScalarGrad  grad_mu = variable_list.template get_gradient<ScalarGrad>(1);
  ScalarValue n1      = variable_list.template get_value<ScalarValue>(2);
  ScalarGrad  grad_n1 = variable_list.template get_gradient<ScalarGrad>(2);
  VectorGrad  grad_u  = variable_list.template get_gradient<VectorGrad>(3);

  ScalarValue h1V  = 3.0 * n1 * n1 - 2.0 * n1 * n1 * n1;
  ScalarValue hn1V = 6.0 * n1 - 6.0 * n1 * n1;

  ScalarValue c_alpha =
    ((B2 * c + 0.5 * (B1 - A1) * h1V) / (A2 * h1V + B2 * (1.0 - h1V)));
  ScalarValue c_beta =
    ((A2 * c + 0.5 * (A1 - B1) * (1.0 - h1V)) / (A2 * h1V + B2 * (1.0 - h1V)));

  ScalarValue faV   = (A2 * c_alpha * c_alpha + A1 * c_alpha + A0);
  ScalarValue facV  = (2.0 * A2 * c_alpha + A1);
  ScalarValue faccV = (2.0 * A2);
  ScalarValue fbV   = (B2 * c_beta * c_beta + B1 * c_beta + B0);
  ScalarValue fbccV = (2.0 * B2);

  ScalarValue fbarriernV = (2.0 * n1 - 6.0 * n1 * n1 + 4.0 * n1 * n1 * n1);

  ScalarValue cacV, canV, cbnV, cbcV, cbcnV;

  cacV = fbccV / ((1.0 - h1V) * fbccV + h1V * faccV);
  canV = hn1V * (c_alpha - c_beta) * cacV;

  cbcV  = faccV / ((1.0 - h1V) * fbccV + h1V * faccV);
  cbnV  = hn1V * (c_alpha - c_beta) * cbcV;
  cbcnV = (faccV * (fbccV - faccV) * hn1V) /
          (((1.0 - h1V) * fbccV + h1V * faccV) * ((1.0 - h1V) * fbccV + h1V * faccV));
  VectorGrad sfts1, sfts1c, sfts1cc, sfts1n, sfts1cn;

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c_beta + b_p
          sfts1[i][j]   = sfts_linear1[i][j] * c_beta + sfts_const1[i][j];
          sfts1c[i][j]  = sfts_linear1[i][j] * cbcV;
          sfts1cc[i][j] = 0.0;
          sfts1n[i][j]  = sfts_linear1[i][j] * cbnV;
          sfts1cn[i][j] = sfts_linear1[i][j] * cbcnV;
        }
    }

  // compute E2=(E-E0)
  ScalarValue E2[dim][dim], S[dim][dim];

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E2[i][j] = 0.5 * (ux[i][j] + ux[j][i]) - (sfts1[i][j] * h1V);
        }
    }

  // compute stress
  // S=C*(E-E0)
  //  Compute stress tensor (which is equal to the residual, Rux)
  ScalarValue CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

  if (n_dependent_stiffness == true)
    {
      for (unsigned int i = 0; i < 2 * dim - 1 + dim / 3; i++)
        {
          for (unsigned int j = 0; j < 2 * dim - 1 + dim / 3; j++)
            {
              CIJ_combined[i][j] = CIJ_Mg[i][j] * (1.0 - h1V) + CIJ_Beta[i][j] * h1V;
            }
        }
      computeStress<dim>(CIJ_combined, E2, S);
    }
  else
    {
      computeStress<dim>(CIJ_Mg, E2, S);
    }

  // Compute one of the stress terms in the order parameter chemical potential,
  // nDependentMisfitACp = -C*(E-E0)*(E0_n)
  ScalarValue nDependentMisfitAC1 = 0.0;

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          nDependentMisfitAC1 += -S[i][j] * (sfts1n[i][j] * h1V + sfts1[i][j] * hn1V);
        }
    }

  // Compute the other stress term in the order parameter chemical potential,
  // heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
  ScalarValue heterMechAC1 = 0.0;
  ScalarValue S2[dim][dim];

  if (n_dependent_stiffness == true)
    {
      computeStress<dim>(CIJ_Beta - CIJ_Mg, E2, S2);

      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              heterMechAC1 += S2[i][j] * E2[i][j];
            }
        }
      heterMechAC1 = 0.5 * hn1V * heterMechAC1;
    }

  // compute K*nx
  ScalarGrad Knx1;
  for (unsigned int a = 0; a < dim; a++)
    {
      Knx1[a] = 0.0;
      for (unsigned int b = 0; b < dim; b++)
        {
          Knx1[a] += Kn1[a][b] * n1x[b];
        }
    }

  // The terms in the governing equations
  ScalarValue eq_c  = (c);
  ScalarGrad  eqx_c = (-userInputs.dtValue * McV * (h1V * faccV + (1.0 - h1V) * fbccV) /
                      (faccV * fbccV) * mux);

  ScalarValue eq_n1  = (n1 - userInputs.dtValue * Mn1V *
                              ((fbV - faV) * hn1V - (c_beta - c_alpha) * facV * hn1V +
                               W * fbarriernV + nDependentMisfitAC1 + heterMechAC1));
  ScalarGrad  eqx_n1 = (-userInputs.dtValue * Mn1V * Knx1);

  // --- Submitting the terms for the governing equations ---

  variable_list.set_value_term(0, eq_c);
  variable_list.set_gradient_term(0, eqx_c);

  variable_list.set_value_term(2, eq_n1);
  variable_list.set_gradient_term(2, eqx_n1);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           index) const
{
  // --- Getting the values and derivatives of the model variables ---

  // The concentration and its derivatives (names here should match those in the
  // macros above)
  ScalarValue c = variable_list.template get_value<ScalarValue>(0);

  // The first order parameter and its derivatives (names here should match
  // those in the macros above)
  ScalarValue n1 = variable_list.template get_value<ScalarValue>(2);

  // The derivative of the displacement vector (names here should match those in
  // the macros above)
  VectorGrad ux = variable_list.template get_gradient<VectorGrad>(3);

  // --- Setting the expressions for the terms in the governing equations ---

  ScalarValue h1V  = (3.0 * n1 * n1 - 2.0 * n1 * n1 * n1);
  ScalarValue hn1V = (6.0 * n1 - 6.0 * n1 * n1);

  // Calculate c_alpha and c_beta from c
  ScalarValue c_alpha =
    ((B2 * c + 0.5 * (B1 - A1) * h1V) / (A2 * h1V + B2 * (1.0 - h1V)));
  ScalarValue c_beta =
    ((A2 * c + 0.5 * (A1 - B1) * (1.0 - h1V)) / (A2 * h1V + B2 * (1.0 - h1V)));

  ScalarValue facV  = (2.0 * A2 * c_alpha + A1);
  ScalarValue faccV = (constV(2.0) * A2);
  ScalarValue fbcV  = (2.0 * B2 * c_beta + B1);
  ScalarValue fbccV = (constV(2.0) * B2);

  // This double-well function can be used to tune the interfacial energy
  // ScalarValue fbarrierV = (n1*n1-2.0*n1*n1*n1+n1*n1*n1*n1);
  // ScalarValue fbarriernV = (2.0*n1-6.0*n1*n1+4.0*n1*n1*n1);

  // Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
  ScalarValue cacV, canV, cbnV, cbcV, cbcnV;

  cacV = fbccV / ((constV(1.0) - h1V) * fbccV + h1V * faccV);
  canV = hn1V * (c_alpha - c_beta) * cacV;

  cbcV  = faccV / ((constV(1.0) - h1V) * fbccV + h1V * faccV);
  cbnV  = hn1V * (c_alpha - c_beta) * cbcV;
  cbcnV = (faccV * (fbccV - faccV) * hn1V) /
          (((1.0 - h1V) * fbccV + h1V * faccV) *
           ((1.0 - h1V) * fbccV +
            h1V * faccV)); // Note: this is only true if faV and fbV are quadratic

  // Calculate the stress-free transformation strain and its derivatives at the
  // quadrature point
  Tensor<2, dim, VectorizedArray<double>> sfts1, sfts1c, sfts1cc, sfts1n, sfts1cn;

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          // Polynomial fits for the stress-free transformation strains, of the
          // form: sfts = a_p * c_beta + b_p
          sfts1[i][j]   = constV(sfts_linear1[i][j]) * c_beta + constV(sfts_const1[i][j]);
          sfts1c[i][j]  = constV(sfts_linear1[i][j]) * cbcV;
          sfts1cc[i][j] = constV(0.0);
          sfts1n[i][j]  = constV(sfts_linear1[i][j]) * cbnV;
          sfts1cn[i][j] = constV(sfts_linear1[i][j]) * cbcnV;
        }
    }

  // compute E2=(E-E0)
  VectorizedArray<double> E2[dim][dim], S[dim][dim];

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          E2[i][j] = constV(0.5) * (ux[i][j] + ux[j][i]) - (sfts1[i][j] * h1V);
        }
    }

  // compute stress
  // S=C*(E-E0)
  //  Compute stress tensor (which is equal to the residual, Rux)
  VectorizedArray<double> CIJ_combined[CIJ_tensor_size][CIJ_tensor_size];

  if (n_dependent_stiffness == true)
    {
      for (unsigned int i = 0; i < 2 * dim - 1 + dim / 3; i++)
        {
          for (unsigned int j = 0; j < 2 * dim - 1 + dim / 3; j++)
            {
              CIJ_combined[i][j] =
                CIJ_Mg[i][j] * (constV(1.0) - h1V) + CIJ_Beta[i][j] * h1V;
            }
        }
      computeStress<dim>(CIJ_combined, E2, S);
    }
  else
    {
      computeStress<dim>(CIJ_Mg, E2, S);
    }

  VectorGrad eqx_u;

  // Fill residual corresponding to mechanics
  // R=-C*(E-E0)

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          eqx_u[i][j] = -S[i][j];
        }
    }

  ScalarValue mu_c = constV(0.0);
  mu_c += facV * cacV * (constV(1.0) - h1V) + fbcV * cbcV * h1V;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          mu_c -= S[i][j] * (sfts1c[i][j] * h1V);
        }
    }

  ScalarValue eq_mu = (mu_c);

  // --- Submitting the terms for the governing equations ---

  variable_list.set_scalar_value_term_rhs(1, eq_mu);
  variable_list.set_vector_gradient_term_rhs(3, eqx_u);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           index) const
{
  // --- Getting the values and derivatives of the model variables ---

  // n1
  ScalarValue n1 = variable_list.template get_value<ScalarValue>(2);

  // --- Setting the expressions for the terms in the governing equations ---

  VectorGrad eqx_Du;

  // Take advantage of E being simply 0.5*(ux + transpose(ux)) and use the
  // dealii "symmetrize" function
  Tensor<2, dim, VectorizedArray<double>> E;
  E = symmetrize(variable_list.get_change_in_vector_gradient(3));

  // Compute stress tensor (which is equal to the residual, Rux)
  if (n_dependent_stiffness == true)
    {
      ScalarValue h1V = (3.0 * n1 * n1 - 2.0 * n1 * n1 * n1);
      Tensor<2, CIJ_tensor_size, VectorizedArray<double>> CIJ_combined;
      CIJ_combined = CIJ_Mg * (constV(1.0) - h1V);
      CIJ_combined += CIJ_Beta * (h1V);

      computeStress<dim>(CIJ_combined, E, eqx_Du);
    }
  else
    {
      computeStress<dim>(CIJ_Mg, E, eqx_Du);
    }

  // --- Submitting the terms for the governing equations ---

  variable_list.set_vector_gradient_term_lhs(3, eqx_Du);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
