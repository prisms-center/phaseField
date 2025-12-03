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
  set_variable_type(0, FieldInfo::TensorRank::Vector);
  set_variable_equation_type(0, ImplicitTimeDependent);

  set_dependencies_value_term_rhs(0, "u, u_star, grad(p)");
  set_dependencies_gradient_term_rhs(0, "u, u_star, grad(p), old_1(u)");
  set_dependencies_value_term_lhs(0, "change(u)");
  set_dependencies_gradient_term_lhs(0, "change(u), old_1(u)");
  set_solve_block(0, 1);

  set_variable_name(1, "u_star");
  set_variable_type(1, FieldInfo::TensorRank::Vector);
  set_variable_equation_type(1, ImplicitTimeDependent);

  set_dependencies_value_term_rhs(1, "u_star, u, grad(u), div(u)");
  set_dependencies_gradient_term_rhs(1, "u_star, u, grad(u), lap(u)");
  set_dependencies_value_term_lhs(1, "change(u_star)");
  set_dependencies_gradient_term_lhs(1, "change(u_star), u");
  set_solve_block(1, 0);

  set_variable_name(2, "p");
  set_variable_type(2, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(2, TimeIndependent);

  set_dependencies_value_term_rhs(2, "div(u_star)");
  set_dependencies_gradient_term_rhs(2, "grad(p), u_star, u, grad(u), lap(u)");
  set_dependencies_value_term_lhs(2, "");
  set_dependencies_gradient_term_lhs(2, "grad(change(p))");
  set_solve_block(2, 0);
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
  if (index == 1 && solve_block == 0)
    {
      VectorValue u      = variable_list.template get_value<VectorValue>(0);
      VectorGrad  grad_u = variable_list.template get_gradient<VectorGrad>(0);
      ScalarValue div_u  = variable_list.template get_divergence<ScalarValue>(0);
      VectorValue lap_u  = variable_list.template get_laplacian<VectorValue>(0);
      VectorValue u_star = variable_list.template get_value<VectorValue>(1);

      ScalarValue stabilization_parameter =
        compute_stabilization_parameter(u, element_volume);

      VectorValue advection_term;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              advection_term[i] += u[j] * grad_u[i][j];
            }
        }
      // Add the skew-symmetric contribution
      advection_term += 0.5 * div_u * u;

      // Calculate the residual
      ScalarGrad residual = (u - u_star) - (dt * advection_term) + (dt * nu * lap_u);

      // Compute the SUPG term which is the dot product of the the velocity and the
      // gradient of the shape function. The shape function is not directly accessible in
      // PRISMS-PF, so we create a dummy VectorGrad and fill in the diagonal to achieve
      // the same result.
      VectorGrad supg_term;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              supg_term[i][j] = residual[i] * u[j];
            }
        }
      supg_term *= stabilization_parameter;

      VectorValue eq_u      = u - u_star - (dt * advection_term);
      VectorGrad  eq_grad_u = (-dt * grad_u * nu) + supg_term;

      variable_list.set_value_term(1, eq_u);
      variable_list.set_gradient_term(1, eq_grad_u);
    }
  if (index == 2 && solve_block == 0)
    {
      VectorValue u          = variable_list.template get_value<VectorValue>(0);
      VectorGrad  grad_u     = variable_list.template get_gradient<VectorGrad>(0);
      VectorValue lap_u      = variable_list.template get_laplacian<VectorValue>(0);
      VectorValue u_star     = variable_list.template get_value<VectorValue>(1);
      ScalarValue div_u_star = variable_list.template get_divergence<ScalarValue>(1);
      ScalarGrad  grad_p     = variable_list.template get_gradient<ScalarGrad>(2);

      ScalarValue stabilization_parameter =
        compute_stabilization_parameter(u_star, element_volume);

      ScalarValue eq_p      = -div_u_star / (dt + stabilization_parameter);
      ScalarGrad  eq_grad_p = -grad_p;

      VectorValue advection_term;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              advection_term[i] += u[j] * grad_u[i][j];
            }
        }
      // Add the skew-symmetric contribution
      advection_term += 0.5 * div_u_star * u_star;

      // Calculate the residual
      ScalarGrad residual = ((u_star - u) / dt) + advection_term - (nu * lap_u);

      // Add the residual to the poisson solve
      eq_grad_p -= residual * stabilization_parameter / (dt + stabilization_parameter);

      variable_list.set_value_term(2, eq_p);
      variable_list.set_gradient_term(2, eq_grad_p);
    }
  if (index == 0 && solve_block == 1)
    {
      VectorValue u       = variable_list.template get_value<VectorValue>(0);
      VectorValue u_old_1 = variable_list.template get_value<VectorValue>(0, OldOne);
      VectorValue u_star  = variable_list.template get_value<VectorValue>(1);
      ScalarGrad  grad_p  = variable_list.template get_gradient<ScalarGrad>(2);

      ScalarValue stabilization_parameter =
        compute_stabilization_parameter(u_old_1, element_volume);

      // Calculate the residual
      ScalarGrad residual = u_star - u - (dt * grad_p);

      // Compute the SUPG term which is the dot product of the the velocity and the
      // gradient of the shape function. The shape function is not directly accessible in
      // PRISMS-PF, so we create a dummy VectorGrad and fill in the diagonal to achieve
      // the same result.
      VectorGrad eq_grad_u;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              eq_grad_u[i][j] = residual[i] * u_old_1[j];
            }
        }
      eq_grad_u *= stabilization_parameter;

      VectorValue eq_u = u_star - u - (dt * grad_p);

      variable_list.set_value_term(0, eq_u);
      variable_list.set_gradient_term(0, eq_grad_u);
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
  if (index == 1 && solve_block == 0)
    {
      VectorValue u = variable_list.template get_value<VectorValue>(0);
      VectorValue change_u_star =
        variable_list.template get_value<VectorValue>(1, Change);

      ScalarValue stabilization_parameter =
        compute_stabilization_parameter(u, element_volume);

      // Compute the SUPG term which is the dot product of the the velocity and the
      // gradient of the shape function. The shape function is not directly accessible in
      // PRISMS-PF, so we create a dummy VectorGrad and fill in the diagonal to achieve
      // the same result.
      VectorGrad eq_change_grad_u_star;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              eq_change_grad_u_star[i][j] = change_u_star[i] * u[j];
            }
        }
      eq_change_grad_u_star *= stabilization_parameter;

      variable_list.set_value_term(1, change_u_star, Change);
      variable_list.set_gradient_term(1, eq_change_grad_u_star, Change);
    }
  if (index == 2 && solve_block == 0)
    {
      ScalarGrad change_grad_p =
        variable_list.template get_gradient<ScalarGrad>(2, Change);

      variable_list.set_gradient_term(2, change_grad_p, Change);
    }
  if (index == 0 && solve_block == 1)
    {
      VectorValue u_old_1  = variable_list.template get_value<VectorValue>(0, OldOne);
      VectorValue change_u = variable_list.template get_value<VectorValue>(0, Change);

      ScalarValue stabilization_parameter =
        compute_stabilization_parameter(u_old_1, element_volume);

      // Compute the SUPG term which is the dot product of the the velocity and the
      // gradient of the shape function. The shape function is not directly accessible in
      // PRISMS-PF, so we create a dummy VectorGrad and fill in the diagonal to achieve
      // the same result.
      VectorGrad eq_change_grad_u;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              eq_change_grad_u[i][j] = change_u[i] * u_old_1[j];
            }
        }
      eq_change_grad_u *= stabilization_parameter;

      variable_list.set_value_term(0, change_u, Change);
      variable_list.set_gradient_term(0, eq_change_grad_u, Change);
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
