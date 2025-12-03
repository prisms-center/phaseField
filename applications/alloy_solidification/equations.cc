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
  set_variable_name(0, "U");
  set_variable_type(0, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(0, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(0, "U,xi,phi,grad(phi),grad(U)");
  set_dependencies_gradient_term_rhs(0, "U,grad(U),grad(phi),phi,xi");

  set_variable_name(1, "phi");
  set_variable_type(1, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(1, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(1, "phi,U,xi");
  set_dependencies_gradient_term_rhs(1, "");

  set_variable_name(2, "xi");
  set_variable_type(2, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(2, Auxiliary);
  set_dependencies_value_term_rhs(2, "phi,U,grad(phi)");
  set_dependencies_gradient_term_rhs(2, "grad(phi)");

  set_variable_name(3, "c");
  set_variable_type(3, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(3, ExplicitTimeDependent);
  set_is_postprocessed_field(3, true);
  set_dependencies_value_term_rhs(3, "phi,U");
  set_dependencies_gradient_term_rhs(3, "");
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  // --- Getting the values and derivatives of the model variables ---

  // The dimensionless solute supersaturation and its derivatives
  ScalarValue U  = variable_list.template get_value<ScalarValue>(0);
  ScalarGrad  Ux = variable_list.template get_gradient<ScalarGrad>(0);

  // The order parameter and its derivatives
  ScalarValue phi  = variable_list.template get_value<ScalarValue>(1);
  ScalarGrad  phix = variable_list.template get_gradient<ScalarGrad>(1);

  // The auxiliary parameter and its derivatives
  ScalarValue xi = variable_list.template get_value<ScalarValue>(2);

  // --- Setting the expressions for the terms in the governing equations ---

  // Calculation of interface normal vector
  ScalarValue normgradn = std::sqrt(phix.norm_square());
  ScalarGrad  normal    = phix / (normgradn + regval);

  // The cosine of theta
  ScalarValue cth = normal[0];
  // The sine of theta
  ScalarValue sth = normal[1];
  // The cosine of 4 theta
  ScalarValue c4th =
    (sth * sth * sth * sth) + (cth * cth * cth * cth) - (6.0 * sth * sth * cth * cth);

  // Anisotropic term
  ScalarValue a_n = 1.0 + (epsilon * c4th);

  // coeffcient before phi
  ScalarValue tau_phi = (1.0 + (1.0 - k) * U) * a_n * a_n;

  // coeffcient before U
  ScalarValue tau_U = (((1.0 + k) / 2.0) - ((1.0 - k) * phi / 2.0));

  // Antitrapping term
  ScalarGrad j_at;
  j_at = (1.0 / (2.0 * std::sqrt(2.0))) * (1.0 + (1.0 - k) * U) * (xi / tau_phi) * normal;

  // grad_phi and grad_U dot product term
  ScalarValue val_term1 =
    get_timestep() * (1.0 + (1.0 - k) * U) * xi / (2.0 * tau_phi * tau_U);
  ScalarValue val_term2 = get_timestep() * ((1.0 - k) / 2.0) / (tau_U * tau_U) *
                          (phix[0] * (Dtilde * ((1.0 - phi) / 2.0) * Ux[0] + j_at[0]) +
                           phix[1] * (Dtilde * ((1.0 - phi) / 2.0) * Ux[1] + j_at[1]));

  // Define required equations
  ScalarValue eq_U = (U + val_term1 - val_term2);

  ScalarGrad eqx_U =
    (-1.0 * get_timestep() * (Dtilde * ((1.0 - phi) / 2.0) * Ux + j_at) / tau_U);

  ScalarValue eq_phi = phi + (get_timestep() * xi / tau_phi);

  // --- Submitting the terms for the governing equations ---

  // Terms for the equation to evolve the concentration
  variable_list.set_value_term(0, eq_U);
  variable_list.set_gradient_term(0, eqx_U);

  // Terms for the equation to evolve the order parameter
  variable_list.set_value_term(1, eq_phi);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           current_index) const
{
  if (current_index == 2)
    {
      // --- Getting the values and derivatives of the model variables ---

      ScalarValue U    = variable_list.template get_value<ScalarValue>(0);
      ScalarValue phi  = variable_list.template get_value<ScalarValue>(1);
      ScalarGrad  phix = variable_list.template get_gradient<ScalarGrad>(1);

      // --- Setting the expressions for the terms in the governing equations ---

      // Calculation of interface normal vector
      ScalarValue normgradn = std::sqrt(phix.norm_square());
      ScalarGrad  normal    = phix / (normgradn + regval);

      // The cosine of theta
      ScalarValue cth = normal[0];
      // The sine of theta
      ScalarValue sth = normal[1];

      // The cosine of 4 theta
      ScalarValue c4th =
        (sth * sth * sth * sth) + (cth * cth * cth * cth) - (6.0 * sth * sth * cth * cth);
      // The sine of 4 theta
      ScalarValue s4th = (4.0 * sth * cth * cth * cth) - (4.0 * sth * sth * sth * cth);

      // Anisotropic term
      ScalarValue a_n = 1.0 + (epsilon * c4th);

      // gradient energy coefficient, its derivative and square
      ScalarValue a_d = -4.0 * epsilon * s4th;

      // dimensionless temperature changes
      ScalarValue y = q_point_loc[1]; // The y-component
      ScalarValue t_n =
        get_user_inputs().get_temporal_discretization().get_time(); // The time
      ScalarValue tep = ((y - y0 - Vtilde * t_n) / ltilde);

      // The anisotropy term that enters in to the equation for xi
      ScalarGrad aniso;
      aniso[0] = a_n * a_n * phix[0] - a_n * a_d * phix[1];
      aniso[1] = a_n * a_n * phix[1] + a_n * a_d * phix[0];

      // Define the terms in the equations
      ScalarValue eq_xi =
        phi - (phi * phi * phi) -
        (lamda * (1.0 - phi * phi) * (1.0 - phi * phi) * (U + tep + U_off));

      ScalarGrad eqx_xi = -aniso;

      // --- Submitting the terms for the governing equations ---

      variable_list.set_value_term(2, eq_xi);
      variable_list.set_gradient_term(2, eqx_xi);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           current_index) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  ScalarValue U   = variable_list.template get_value<ScalarValue>(0);
  ScalarValue phi = variable_list.template get_value<ScalarValue>(1);

  ScalarValue c = (c0 / 2.0 / (1 + U0 - U0 * k)) * ((1.0 + k) - (1.0 - k) * phi) *
                  ((1.0) + (1.0 - k) * U);

  variable_list.set_value_term(3, c);
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE