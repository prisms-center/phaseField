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
  set_variable_type(0, Scalar);
  set_variable_equation_type(0, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(0, "n, n_2, n_3, n_4");
  set_dependencies_gradient_term_rhs(0, "grad(n), grad(n_2), grad(n_3), grad(n_4)");
  set_solve_block(0, 1);

  set_variable_name(1, "n_2");
  set_variable_type(1, Scalar);
  set_variable_equation_type(1, Auxiliary);
  set_dependencies_value_term_rhs(1, "n");
  set_dependencies_gradient_term_rhs(1, "grad(n)");

  set_variable_name(2, "n_3");
  set_variable_type(2, Scalar);
  set_variable_equation_type(2, Auxiliary);
  set_dependencies_value_term_rhs(2, "n, n_2");
  set_dependencies_gradient_term_rhs(2, "grad(n_2)");

  set_variable_name(3, "n_4");
  set_variable_type(3, Scalar);
  set_variable_equation_type(3, Auxiliary);
  set_dependencies_value_term_rhs(3, "n, n_3");
  set_dependencies_gradient_term_rhs(3, "grad(n_3)");
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
      ScalarValue n        = variable_list.template get_value<ScalarValue>(0);
      ScalarGrad  grad_n   = variable_list.template get_gradient<ScalarGrad>(0);
      ScalarValue n_2      = variable_list.template get_value<ScalarValue>(1);
      ScalarGrad  grad_n_2 = variable_list.template get_gradient<ScalarGrad>(1);
      ScalarValue n_3      = variable_list.template get_value<ScalarValue>(2);
      ScalarGrad  grad_n_3 = variable_list.template get_gradient<ScalarGrad>(2);
      ScalarValue n_4      = variable_list.template get_value<ScalarValue>(3);
      ScalarGrad  grad_n_4 = variable_list.template get_gradient<ScalarGrad>(3);

      ScalarValue k_1_bulk_energy = 4.0 * n * (n - 1.0) * (n - 0.5);
      ScalarValue k_2_bulk_energy = 4.0 * n_2 * (n_2 - 1.0) * (n_2 - 0.5);
      ScalarValue k_3_bulk_energy = 4.0 * n_3 * (n_3 - 1.0) * (n_3 - 0.5);
      ScalarValue k_4_bulk_energy = 4.0 * n_4 * (n_4 - 1.0) * (n_4 - 0.5);

      ScalarValue eq_n      = n - (this->get_timestep() * MnV / 6.0 *
                              (k_1_bulk_energy + 2.0 * k_2_bulk_energy +
                               2.0 * k_3_bulk_energy + k_4_bulk_energy));
      ScalarGrad  eq_grad_n = -this->get_timestep() * KnV * MnV / 6.0 *
                             (grad_n + 2.0 * grad_n_2 + 2.0 * grad_n_3 + grad_n_4);

      variable_list.set_value_term(0, eq_n);
      variable_list.set_gradient_term(0, eq_grad_n);
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
  if (solve_block == 0 && index == 1)
    {
      ScalarValue n      = variable_list.template get_value<ScalarValue>(0);
      ScalarGrad  grad_n = variable_list.template get_gradient<ScalarGrad>(0);

      ScalarValue bulk_energy = 4.0 * n * (n - 1.0) * (n - 0.5);
      ScalarValue eq_n_2      = n - 0.5 * this->get_timestep() * MnV * bulk_energy;
      ScalarGrad  eq_grad_n_2 = -0.5 * this->get_timestep() * MnV * KnV * grad_n;

      variable_list.set_value_term(1, eq_n_2);
      variable_list.set_gradient_term(1, eq_grad_n_2);
    }
  if (solve_block == 0 && index == 2)
    {
      ScalarValue n        = variable_list.template get_value<ScalarValue>(0);
      ScalarValue n_2      = variable_list.template get_value<ScalarValue>(1);
      ScalarGrad  grad_n_2 = variable_list.template get_gradient<ScalarGrad>(1);

      ScalarValue bulk_energy = 4.0 * n_2 * (n_2 - 1.0) * (n_2 - 0.5);
      ScalarValue eq_n_3      = n - 0.5 * this->get_timestep() * MnV * bulk_energy;
      ScalarGrad  eq_grad_n_3 = -0.5 * this->get_timestep() * MnV * KnV * grad_n_2;

      variable_list.set_value_term(2, eq_n_3);
      variable_list.set_gradient_term(2, eq_grad_n_3);
    }
  if (solve_block == 0 && index == 3)
    {
      ScalarValue n        = variable_list.template get_value<ScalarValue>(0);
      ScalarValue n_3      = variable_list.template get_value<ScalarValue>(2);
      ScalarGrad  grad_n_3 = variable_list.template get_gradient<ScalarGrad>(2);

      ScalarValue bulk_energy = 4.0 * n_3 * (n_3 - 1.0) * (n_3 - 0.5);
      ScalarValue eq_n_3      = n - this->get_timestep() * MnV * bulk_energy;
      ScalarGrad  eq_grad_n_3 = -this->get_timestep() * MnV * KnV * grad_n_3;

      variable_list.set_value_term(3, eq_n_3);
      variable_list.set_gradient_term(3, eq_grad_n_3);
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
