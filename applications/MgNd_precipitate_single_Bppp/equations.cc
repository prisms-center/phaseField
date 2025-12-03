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
  set_variable_type(0, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(0, ExplicitTimeDependent);

  set_dependencies_value_term_rhs(0, "c");
  set_dependencies_gradient_term_rhs(0, "n1, grad(mu)");

  set_variable_name(1, "mu");
  set_variable_type(1, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(1, Auxiliary);

  set_dependencies_value_term_rhs(1, "c, n1, grad(u)");
  set_dependencies_gradient_term_rhs(1, "");

  set_variable_name(2, "n1");
  set_variable_type(2, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(2, ExplicitTimeDependent);

  set_dependencies_value_term_rhs(2, "c, n1, grad(u)");
  set_dependencies_gradient_term_rhs(2, "grad(n1)");

  set_variable_name(3, "u");
  set_variable_type(3, FieldInfo::TensorRank::Vector);
  set_variable_equation_type(3, TimeIndependent);

  set_dependencies_value_term_rhs(3, "");
  set_dependencies_gradient_term_rhs(3, "c, n1, grad(u)");
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
  VectorGrad  grad_u  = variable_list.template get_symmetric_gradient<VectorGrad>(3);

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

  // Compute strain
  VectorGrad strain = grad_u - (sfts1 * h1V);

  // Compute stress
  VectorGrad stress;
  if (n_dependent_stiffness == true)
    {
      dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> CIJ_combined =
        CIJ_Mg * (1.0 - h1V) + CIJ_Beta * h1V;
      compute_stress<dim, ScalarValue>(CIJ_combined, strain, stress);
    }
  else
    {
      compute_stress<dim, ScalarValue>(CIJ_Mg, strain, stress);
    }

  // Compute one of the stress terms in the order parameter chemical potential,
  ScalarValue nDependentMisfitAC1 = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          nDependentMisfitAC1 +=
            -stress[i][j] * (sfts1n[i][j] * h1V + sfts1[i][j] * hn1V);
        }
    }

  // Compute the other stress term in the order parameter chemical potential,
  ScalarValue heterMechAC1 = 0.0;
  VectorGrad  stress_2;
  if (n_dependent_stiffness == true)
    {
      compute_stress<dim, ScalarValue>(CIJ_Beta - CIJ_Mg, strain, stress_2);

      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              heterMechAC1 += stress_2[i][j] * strain[i][j];
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
          Knx1[a] += Kn1[a][b] * grad_n1[b];
        }
    }

  ScalarValue eq_c  = (c);
  ScalarGrad  eqx_c = (-get_timestep() * McV * (h1V * faccV + (1.0 - h1V) * fbccV) /
                      (faccV * fbccV) * grad_mu);

  ScalarValue eq_n1  = (n1 - get_timestep() * Mn1V *
                              ((fbV - faV) * hn1V - (c_beta - c_alpha) * facV * hn1V +
                               W * fbarriernV + nDependentMisfitAC1 + heterMechAC1));
  ScalarGrad  eqx_n1 = (-get_timestep() * Mn1V * Knx1);

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
  if (index == 1)
    {
      ScalarValue c      = variable_list.template get_value<ScalarValue>(0);
      ScalarValue n1     = variable_list.template get_value<ScalarValue>(2);
      VectorGrad  grad_u = variable_list.template get_symmetric_gradient<VectorGrad>(3);

      ScalarValue h1V  = (3.0 * n1 * n1 - 2.0 * n1 * n1 * n1);
      ScalarValue hn1V = (6.0 * n1 - 6.0 * n1 * n1);

      // Calculate c_alpha and c_beta from c
      ScalarValue c_alpha =
        ((B2 * c + 0.5 * (B1 - A1) * h1V) / (A2 * h1V + B2 * (1.0 - h1V)));
      ScalarValue c_beta =
        ((A2 * c + 0.5 * (A1 - B1) * (1.0 - h1V)) / (A2 * h1V + B2 * (1.0 - h1V)));

      ScalarValue facV  = (2.0 * A2 * c_alpha + A1);
      ScalarValue faccV = (2.0 * A2);
      ScalarValue fbcV  = (2.0 * B2 * c_beta + B1);
      ScalarValue fbccV = (2.0 * B2);

      // Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
      ScalarValue cacV, canV, cbnV, cbcV, cbcnV;

      cacV = fbccV / ((1.0 - h1V) * fbccV + h1V * faccV);
      canV = hn1V * (c_alpha - c_beta) * cacV;

      cbcV  = faccV / ((1.0 - h1V) * fbccV + h1V * faccV);
      cbnV  = hn1V * (c_alpha - c_beta) * cbcV;
      cbcnV = (faccV * (fbccV - faccV) * hn1V) /
              (((1.0 - h1V) * fbccV + h1V * faccV) * ((1.0 - h1V) * fbccV + h1V * faccV));

      // Calculate the stress-free transformation strain and its derivatives at the
      // quadrature point
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

      // Compute strain
      VectorGrad strain = grad_u - (sfts1 * h1V);

      // Compute stress
      VectorGrad stress;
      if (n_dependent_stiffness == true)
        {
          dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> CIJ_combined =
            CIJ_Mg * (1.0 - h1V) + CIJ_Beta * h1V;
          compute_stress<dim, ScalarValue>(CIJ_combined, strain, stress);
        }
      else
        {
          compute_stress<dim, ScalarValue>(CIJ_Mg, strain, stress);
        }

      ScalarValue mu_c = 0.0;
      mu_c += facV * cacV * (1.0 - h1V) + fbcV * cbcV * h1V;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              mu_c -= stress[i][j] * (sfts1c[i][j] * h1V);
            }
        }
      ScalarValue eq_mu = (mu_c);
      variable_list.set_value_term(1, eq_mu);
    }
  if (index == 3)
    {
      ScalarValue c      = variable_list.template get_value<ScalarValue>(0);
      ScalarValue n1     = variable_list.template get_value<ScalarValue>(2);
      VectorGrad  grad_u = variable_list.template get_symmetric_gradient<VectorGrad>(3);

      ScalarValue h1V = (3.0 * n1 * n1 - 2.0 * n1 * n1 * n1);

      // Calculate c_alpha and c_beta from c
      ScalarValue c_beta =
        ((A2 * c + 0.5 * (A1 - B1) * (1.0 - h1V)) / (A2 * h1V + B2 * (1.0 - h1V)));

      // Calculate the stress-free transformation strain and its derivatives at the
      // quadrature point
      VectorGrad sfts1;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              sfts1[i][j] = sfts_linear1[i][j] * c_beta + sfts_const1[i][j];
            }
        }

      // Compute strain
      VectorGrad strain = grad_u - (sfts1 * h1V);

      // Compute stress
      VectorGrad stress;
      if (n_dependent_stiffness == true)
        {
          dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> CIJ_combined =
            CIJ_Mg * (1.0 - h1V) + CIJ_Beta * h1V;
          compute_stress<dim, ScalarValue>(CIJ_combined, strain, stress);
        }
      else
        {
          compute_stress<dim, ScalarValue>(CIJ_Mg, strain, stress);
        }

      variable_list.set_gradient_term(3, -stress);
    }
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
  if (index == 3)
    {
      ScalarValue n1 = variable_list.template get_value<ScalarValue>(2);
      VectorGrad  change_grad_u =
        variable_list.template get_symmetric_gradient<VectorGrad>(3, Change);

      // Compute strain
      VectorGrad strain = change_grad_u;

      // Compute stress
      VectorGrad stress;
      if (n_dependent_stiffness == true)
        {
          ScalarValue h1V = (3.0 * n1 * n1 - 2.0 * n1 * n1 * n1);
          dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> CIJ_combined =
            CIJ_Mg * (1.0 - h1V) + CIJ_Beta * h1V;
          compute_stress<dim, ScalarValue>(CIJ_combined, strain, stress);
        }
      else
        {
          compute_stress<dim, ScalarValue>(CIJ_Mg, strain, stress);
        }

      // --- Submitting the terms for the governing equations ---

      variable_list.set_gradient_term(3, stress, Change);
    }
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
