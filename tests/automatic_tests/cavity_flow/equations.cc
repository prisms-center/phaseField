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
  set_variable_name(0, "u");
  set_variable_type(0, Vector);
  set_variable_equation_type(0, ExplicitTimeDependent);

  set_dependencies_value_term_rhs(0, "u, grad(P)");
  set_dependencies_gradient_term_rhs(0, "");
  set_solve_block(0, 1);

  set_variable_name(1, "u_star");
  set_variable_type(1, Vector);
  set_variable_equation_type(1, ExplicitTimeDependent);

  set_dependencies_value_term_rhs(1, "u, grad(u), div(u)");
  set_dependencies_gradient_term_rhs(1, "grad(u)");
  set_solve_block(1, 0);

  set_variable_name(2, "P");
  set_variable_type(2, Scalar);
  set_variable_equation_type(2, TimeIndependent);

  set_dependencies_value_term_rhs(2, "div(u_star)");
  set_dependencies_gradient_term_rhs(2, "grad(P)");
  set_dependencies_value_term_lhs(2, "");
  set_dependencies_gradient_term_lhs(2, "grad(change(P))");
  set_solve_block(2, 0);

  set_variable_name(3, "continuity");
  set_variable_type(3, Scalar);
  set_variable_equation_type(3, ExplicitTimeDependent);

  set_dependencies_value_term_rhs(3, "div(u)");
  set_dependencies_gradient_term_rhs(3, "");
  set_solve_block(3, 1);
  set_is_postprocessed_field(3, true);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  if (solve_block == 1)
    {
      VectorValue u  = variable_list.template get_value<VectorValue>(0);
      ScalarGrad  Px = variable_list.template get_gradient<ScalarGrad>(2);

      VectorValue eq_u = u - (this->get_timestep() * Px);

      variable_list.set_value_term(0, eq_u);
    }
  if (solve_block == 0)
    {
      VectorValue u     = variable_list.template get_value<VectorValue>(0);
      VectorGrad  ux    = variable_list.template get_gradient<VectorGrad>(0);
      ScalarValue div_u = variable_list.template get_divergence<ScalarValue>(0);

      VectorValue advection_term;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              advection_term[i] += u[j] * ux[i][j];
            }
        }
      // Add the skew-symmetric contribution
      advection_term += 0.5 * div_u * u;

      VectorValue eq_u  = u - (this->get_timestep() * advection_term);
      VectorGrad  eqx_u = -this->get_timestep() / Re * ux;

      variable_list.set_value_term(1, eq_u);
      variable_list.set_gradient_term(1, eqx_u);
    }
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
  if (index == 2)
    {
      ScalarValue div_u_star = variable_list.template get_divergence<ScalarValue>(1);
      ScalarGrad  Px         = variable_list.template get_gradient<ScalarGrad>(2);

      ScalarValue eq_P  = -div_u_star / this->get_timestep();
      ScalarGrad  eqx_P = -Px;

      variable_list.set_value_term(2, eq_P);
      variable_list.set_gradient_term(2, eqx_P);
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
  if (index == 2)
    {
      ScalarGrad change_Px = variable_list.template get_gradient<ScalarGrad>(2, Change);

      variable_list.set_gradient_term(2, change_Px, Change);
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
  ScalarValue div_u = variable_list.template get_divergence<ScalarValue>(0);

  variable_list.set_value_term(3, div_u);
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
