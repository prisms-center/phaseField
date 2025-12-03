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

  set_dependencies_value_term_rhs(0, "n, dndt");
  set_dependencies_gradient_term_rhs(0, "");

  set_variable_name(1, "u");
  set_variable_type(1, FieldInfo::TensorRank::Vector);
  set_variable_equation_type(1, TimeIndependent);

  set_dependencies_value_term_rhs(1, "");
  set_dependencies_gradient_term_rhs(1, "n, grad(u), Ex");
  set_dependencies_value_term_lhs(1, "");
  set_dependencies_gradient_term_lhs(1, "n, grad(change(u)), Ex");

  set_variable_name(2, "dndt");
  set_variable_type(2, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(2, Auxiliary);

  set_dependencies_value_term_rhs(2, "n, grad(u), Ex, Gx");
  set_dependencies_gradient_term_rhs(2, "grad(n), Gx");

  set_variable_name(3, "Ex");
  set_variable_type(3, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(3, Constant);

  set_variable_name(4, "Gx");
  set_variable_type(4, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(4, Constant);

  set_variable_name(5, "f_tot");
  set_variable_type(5, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(5, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(5, "n, grad(n), grad(u), Ex, Gx");
  set_dependencies_gradient_term_rhs(5, "");
  set_is_postprocessed_field(5, true);

  set_variable_name(6, "s11");
  set_variable_type(6, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(6, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(6, "n, grad(u), Ex");
  set_dependencies_gradient_term_rhs(6, "");
  set_is_postprocessed_field(6, true);

  set_variable_name(7, "s12");
  set_variable_type(7, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(7, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(7, "n, grad(u), Ex");
  set_dependencies_gradient_term_rhs(7, "");
  set_is_postprocessed_field(7, true);

  set_variable_name(8, "s22");
  set_variable_type(8, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(8, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(8, "n, grad(u), Ex");
  set_dependencies_gradient_term_rhs(8, "");
  set_is_postprocessed_field(8, true);

  set_variable_name(9, "e22");
  set_variable_type(9, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(9, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(9, "grad(u), Ex");
  set_dependencies_gradient_term_rhs(9, "");
  set_is_postprocessed_field(9, true);

  set_variable_name(10, "f_int");
  set_variable_type(10, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(10, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(10, "n, grad(n), grad(u), Ex, Gx");
  set_dependencies_gradient_term_rhs(10, "");
  set_is_postprocessed_field(10, true);

  set_variable_name(11, "f_el");
  set_variable_type(11, FieldInfo::TensorRank::Scalar);
  set_variable_equation_type(11, ExplicitTimeDependent);
  set_dependencies_value_term_rhs(11, "n, grad(u), Ex");
  set_dependencies_gradient_term_rhs(11, "");
  set_is_postprocessed_field(11, true);
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] const dealii::VectorizedArray<number> &element_volume,
  [[maybe_unused]] Types::Index                           solve_block) const
{
  ScalarValue n    = variable_list.template get_value<ScalarValue>(0);
  ScalarValue dndt = variable_list.template get_value<ScalarValue>(2);

  for (unsigned int j = 0; j < dndt.size(); ++j)
    {
      if (dndt[j] > 0.0)
        {
          dndt[j] = 0.0;
        }
      if (n[j] - dndt[j] * get_timestep() > 1.0)
        {
          dndt[j] = (n[j] - 1.0) / get_timestep();
        }
    }
  ScalarValue eq_n = n - get_timestep() * Mn * dndt;

  variable_list.set_value_term(0, eq_n);
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
      ScalarValue n  = variable_list.template get_value<ScalarValue>(0);
      VectorGrad  ux = variable_list.template get_symmetric_gradient<VectorGrad>(1);
      ScalarValue Ex = variable_list.template get_value<ScalarValue>(3);

      dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> compliance =
        CIJ_base * Ex * (1.0 - 2.0 * n + n * n);
      VectorGrad stress;
      compute_stress<dim, ScalarValue>(compliance, ux, stress);

      variable_list.set_gradient_term(1, -stress);
    }
  if (index == 2)
    {
      ScalarValue n  = variable_list.template get_value<ScalarValue>(0);
      ScalarGrad  nx = variable_list.template get_gradient<ScalarGrad>(0);
      VectorGrad  ux = variable_list.template get_symmetric_gradient<VectorGrad>(1);
      ScalarValue Ex = variable_list.template get_value<ScalarValue>(3);
      ScalarValue Gx = variable_list.template get_value<ScalarValue>(4);

      dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> compliance = CIJ_base * Ex;
      VectorGrad                                             stress;
      compute_stress<dim, ScalarValue>(compliance, ux, stress);
      ScalarValue elastic_energy = 0.0;
      for (unsigned int i = 0; i < dim; i++)
        {
          for (unsigned int j = 0; j < dim; j++)
            {
              elastic_energy += 0.5 * stress[i][j] * ux[i][j];
            }
        }
      ScalarValue eq_n =
        (2.0 * (n - 1.0) * elastic_energy + Gc0 * Gx * 3.0 / 8.0 / ell) * Mn;
      ScalarGrad eqx_n = ell * nx * Gc0 * Gx * 3.0 / 8.0 * Mn;

      variable_list.set_value_term(2, eq_n);
      variable_list.set_gradient_term(2, eqx_n);
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
  if (index == 1)
    {
      ScalarValue n = variable_list.template get_value<ScalarValue>(0);
      VectorGrad  ux_change =
        variable_list.template get_symmetric_gradient<VectorGrad>(1, Change);
      ScalarValue Ex = variable_list.template get_value<ScalarValue>(3);

      dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> compliance =
        CIJ_base * Ex * (1.0 - 2.0 * n + n * n);
      VectorGrad stress;
      compute_stress<dim, ScalarValue>(compliance, ux_change, stress);

      variable_list.set_gradient_term(1, stress, Change);
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
  ScalarValue n  = variable_list.template get_value<ScalarValue>(0);
  ScalarGrad  nx = variable_list.template get_gradient<ScalarGrad>(0);
  VectorGrad  ux = variable_list.template get_symmetric_gradient<VectorGrad>(1);
  ScalarValue Ex = variable_list.template get_value<ScalarValue>(3);
  ScalarValue Gx = variable_list.template get_value<ScalarValue>(4);

  ScalarValue f_int = Gc0 * n * Gx * 3.0 / 8.0 / ell + Gc0 * Gx * 0.5 * ell * nx * nx;

  dealii::Tensor<2, voigt_tensor_size<dim>, ScalarValue> compliance =
    CIJ_base * Ex * (1.0 - 2.0 * n + n * n);
  VectorGrad stress;
  compute_stress<dim, ScalarValue>(compliance, ux, stress);

  ScalarValue f_el = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      for (unsigned int j = 0; j < dim; j++)
        {
          f_el += 0.5 * stress[i][j] * ux[i][j];
        }
    }

  ScalarValue f_tot = f_el + f_int;
  ScalarValue s11   = stress[0][0];
  ScalarValue s12   = stress[0][1];
  ScalarValue s22   = stress[1][1];
  ScalarValue e22   = ux[1][1];

  variable_list.set_value_term(5, f_tot);
  variable_list.set_value_term(6, s11);
  variable_list.set_value_term(7, s12);
  variable_list.set_value_term(8, s22);
  variable_list.set_value_term(9, e22);
  variable_list.set_value_term(10, f_int);
  variable_list.set_value_term(11, f_el);
}

#include "custom_pde.inst"

PRISMS_PF_END_NAMESPACE
