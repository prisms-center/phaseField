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
  set_variable_name(0, "T");
  set_variable_type(0, Scalar);
  set_variable_equation_type(0, TimeIndependent);
  set_dependencies_value_term_rhs(0, "q");
  set_dependencies_gradient_term_rhs(0, "grad(T)");
  set_dependencies_gradient_term_lhs(0, "grad(change(T))");

  set_variable_name(1, "q");
  set_variable_type(1, Scalar);
  set_variable_equation_type(1, Constant);

  set_variable_name(2, "error");
  set_variable_type(2, Scalar);
  set_variable_equation_type(2, ExplicitTimeDependent);
  set_is_postprocessed_field(2, true);
  set_dependencies_value_term_rhs(2, "T");
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index current_index) const
{
  if (current_index == 0)
    {
      ScalarGrad  Tx = variable_list.template get_gradient<Scalar>(0);
      ScalarValue q  = variable_list.template get_value<Scalar>(1);

      variable_list.template set_value_term<Scalar>(0, q);
      variable_list.template set_gradient_term<Scalar>(0, -Tx);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_nonexplicit_lhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc,
  [[maybe_unused]] Types::Index current_index) const
{
  if (current_index == 0)
    {
      ScalarGrad change_Tx = variable_list.template get_gradient<Scalar>(0, Change);

      variable_list.template set_gradient_term<Scalar>(0, change_Tx, Change);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
CustomPDE<dim, degree, number>::compute_postprocess_explicit_rhs(
  [[maybe_unused]] VariableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  ScalarValue T = variable_list.template get_value<Scalar>(0);

  ScalarValue analytic =
    std::sin(M_PI * q_point_loc[0] /
             this->get_user_inputs().get_spatial_discretization().get_size()[0]) *
    q_point_loc[1] / this->get_user_inputs().get_spatial_discretization().get_size()[1] *
    (1.0 -
     q_point_loc[1] / this->get_user_inputs().get_spatial_discretization().get_size()[1]);

  ScalarValue error = (T - analytic) * (T - analytic);

  variable_list.template set_value_term<Scalar>(2, error);
}

INSTANTIATE_TRI_TEMPLATE(CustomPDE)

PRISMS_PF_END_NAMESPACE