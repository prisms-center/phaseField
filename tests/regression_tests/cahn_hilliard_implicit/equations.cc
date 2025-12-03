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
  set_variable_equation_type(0, ImplicitTimeDependent);
  set_dependencies_value_term_rhs(0, "c, old_1(c)");
  set_dependencies_gradient_term_rhs(0, "grad(gamma), c, grad(c)");
  set_dependencies_value_term_lhs(0, "change(c)");
  set_dependencies_gradient_term_lhs(0, "change(c), grad(change(c)), c, grad(c)");

  set_variable_name(1, "gamma");
  set_variable_type(1, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(1, Auxiliary);
  set_dependencies_value_term_rhs(1, "");
  set_dependencies_gradient_term_rhs(1, "grad(c)");

  set_variable_name(2, "f_tot");
  set_variable_type(2, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(2, ExplicitTimeDependent);
  set_is_postprocessed_field(2, true);
  set_dependencies_value_term_rhs(2, "c, grad(c)");
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block,
  [[maybe_unused]] Types::Index                           index) const
{
  if (index == 0)
    {
      ScalarValue c          = variable_list.template get_value<ScalarValue>(0);
      ScalarGrad  grad_c     = variable_list.template get_gradient<ScalarGrad>(0);
      ScalarValue old_c      = variable_list.template get_value<ScalarValue>(0, OldOne);
      ScalarGrad  grad_gamma = variable_list.template get_gradient<ScalarGrad>(1);

      ScalarValue eq_c = old_c - c;
      ScalarGrad  eq_grad_c =
        McV * get_timestep() * (grad_gamma - (12.0 * c * c - 12.0 * c + 2.0) * grad_c);

      variable_list.set_value_term(0, eq_c);
      variable_list.set_gradient_term(0, eq_grad_c);
    }
  if (index == 1)
    {
      ScalarGrad cx = variable_list.template get_gradient<ScalarGrad>(0);

      ScalarGrad eqx_gamma = -KcV * cx;

      variable_list.set_gradient_term(1, eqx_gamma);
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
  if (index == 0)
    {
      ScalarValue change_c = variable_list.template get_value<ScalarValue>(0, Change);
      ScalarGrad  change_grad_c =
        variable_list.template get_gradient<ScalarGrad>(0, Change);
      ScalarValue c      = variable_list.template get_value<ScalarValue>(0);
      ScalarGrad  grad_c = variable_list.template get_gradient<ScalarGrad>(0);

      ScalarValue eq_c = change_c;
      ScalarGrad  eq_grad_c =
        McV * get_timestep() *
        ((12.0 * c * c + 24.0 * c * change_c + 12.0 * change_c * change_c -
          12.0 * (c + change_c) + 2.0) *
           change_grad_c +
         (12.0 * change_c * change_c + 24.0 * c * change_c - 12.0 * change_c) * grad_c);

      variable_list.set_value_term(0, eq_c, Change);
      variable_list.set_gradient_term(0, eq_grad_c, Change);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  ScalarValue c  = variable_list.template get_value<ScalarValue>(0);
  ScalarGrad  cx = variable_list.template get_gradient<ScalarGrad>(0);

  ScalarValue f_tot  = 0.0;
  ScalarValue f_chem = c * c * c * c - 2.0 * c * c * c + c * c;
  ScalarValue f_grad = 0.0;

  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_grad += 0.5 * KcV * cx[i] * cx[j];
        }
    }
  f_tot = f_chem + f_grad;
  variable_list.set_value_term(2, f_tot);
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
