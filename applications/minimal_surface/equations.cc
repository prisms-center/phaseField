#include "custom_pde.h"

#include <prismspf/config.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>

PRISMS_PF_BEGIN_NAMESPACE

void
customAttributeLoader::loadVariableAttributes()
{
  set_variable_name(0, "u");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, TIME_INDEPENDENT);

  set_dependencies_value_term_RHS(0, "");
  set_dependencies_gradient_term_RHS(0, "grad(u)");
  set_dependencies_value_term_LHS(0, "");
  set_dependencies_gradient_term_LHS(0, "grad(change(u)), grad(u)");
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  if (current_index == 0)
    {
      scalarGrad ux = variable_list.get_scalar_gradient(0);

      scalarValue coefficient = 1.0 / std::sqrt(1.0 + ux.norm_square());

      variable_list.set_scalar_gradient_term(0, -coefficient * ux);
    }
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_LHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  if (current_index == 0)
    {
      scalarGrad ux        = variable_list.get_scalar_gradient(0);
      scalarGrad change_ux = variable_list.get_scalar_gradient(0, CHANGE);

      scalarValue coefficient = 1.0 / std::sqrt(1.0 + ux.norm_square());

      scalarGrad term_1 = coefficient * change_ux;
      scalarGrad term_2 = -coefficient * (ux * change_ux) * ux;

      variable_list.set_scalar_gradient_term(0, term_1 + term_2, CHANGE);
    }
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_postprocess_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{}

INSTANTIATE_TRI_TEMPLATE(customPDE)

PRISMS_PF_END_NAMESPACE