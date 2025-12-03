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
  set_dependencies_gradient_term_rhs(0, "n,grad(c),grad(n)");

  set_variable_name(1, "n");
  set_variable_type(1, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(1, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(1, "c,n");
  set_dependencies_gradient_term_rhs(1, "grad(n)");
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

  ScalarValue c  = variable_list.template get_value<ScalarValue>(0);
  ScalarGrad  cx = variable_list.template get_gradient<ScalarGrad>(0);

  ScalarValue n  = variable_list.template get_value<ScalarValue>(1);
  ScalarGrad  nx = variable_list.template get_gradient<ScalarGrad>(1);

  // --- Setting the expressions for the terms in the governing equations ---

  // Free energy for each phase and their first and second derivatives
  ScalarValue fa =
    (-1.6704 - 4.776 * c + 5.1622 * c * c - 2.7375 * c * c * c + 1.3687 * c * c * c * c);
  ScalarValue fac  = (-4.776 + 10.3244 * c - 8.2125 * c * c + 5.4748 * c * c * c);
  ScalarValue facc = (10.3244 - 16.425 * c + 16.4244 * c * c);
  ScalarValue fb   = (5.0 * c * c - 5.9746 * c - 1.5924);
  ScalarValue fbc  = (10.0 * c - 5.9746);
  ScalarValue fbcc = (10.0);

  // Interpolation function and its derivative
  ScalarValue h  = (10.0 * n * n * n - 15.0 * n * n * n * n + 6.0 * n * n * n * n * n);
  ScalarValue hn = (30.0 * n * n - 60.0 * n * n * n + 30.0 * n * n * n * n);

  // Residual equations
  ScalarGrad  mux   = (cx * ((1.0 - h) * facc + h * fbcc) + nx * ((fbc - fac) * hn));
  ScalarValue eq_c  = c;
  ScalarGrad  eqx_c = ((-Mc * get_timestep()) * mux);
  ScalarValue eq_n  = (n - (get_timestep() * Mn) * (fb - fa) * hn);
  ScalarGrad  eqx_n = ((-get_timestep() * Kn * Mn) * nx);

  // --- Submitting the terms for the governing equations ---

  // Terms for the equation to evolve the concentration
  variable_list.set_value_term(0, eq_c);
  variable_list.set_gradient_term(0, eqx_c);

  // Terms for the equation to evolve the order parameter
  variable_list.set_value_term(1, eq_n);
  variable_list.set_gradient_term(1, eqx_n);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           index) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           index) const
{}

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
