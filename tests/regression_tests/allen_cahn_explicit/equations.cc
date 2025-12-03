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
  set_variable_name(0, "n");
  set_variable_type(0, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(0, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(0, "n");
  set_dependencies_gradient_term_rhs(0, "grad(n)");

  set_variable_name(1, "mg_n");
  set_variable_type(1, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(1, ExplicitTimeDependent);
  set_is_postprocessed_field(1, true);
  set_dependencies_value_term_rhs(1, "grad(n)");

  set_variable_name(2, "f_tot");
  set_variable_type(2, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(2, ExplicitTimeDependent);
  set_is_postprocessed_field(2, true);
  set_dependencies_value_term_rhs(2, "n, grad(n)");
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  ScalarValue n  = variable_list.template get_value<ScalarValue>(0);
  ScalarGrad  nx = variable_list.template get_gradient<ScalarGrad>(0);

  ScalarValue fnV   = 4.0 * n * (n - 1.0) * (n - 0.5);
  ScalarValue eq_n  = n - get_timestep() * MnV * fnV;
  ScalarGrad  eqx_n = -get_timestep() * KnV * MnV * nx;

  variable_list.set_value_term(0, eq_n);
  variable_list.set_gradient_term(0, eqx_n);
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
{
  ScalarValue n  = variable_list.template get_value<ScalarValue>(0);
  ScalarGrad  nx = variable_list.template get_gradient<ScalarGrad>(0);

  ScalarValue f_tot  = 0.0;
  ScalarValue f_chem = n * n * n * n - 2.0 * n * n * n + n * n;
  ScalarValue f_grad = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_grad += 0.5 * KnV * nx[i] * nx[j];
        }
    }
  f_tot = f_chem + f_grad;

  variable_list.set_value_term(1, std::sqrt(nx[0] * nx[0] + nx[1] * nx[1]));
  variable_list.set_value_term(2, f_tot);
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
