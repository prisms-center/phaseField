#include "custom_pde.h"

#include <prismspf/config.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attribute_loader.h>

PRISMS_PF_BEGIN_NAMESPACE

void
customAttributeLoader::loadVariableAttributes()
{
  set_variable_name(0, "T");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, TIME_INDEPENDENT);
  set_dependencies_value_term_RHS(0, "q");
  set_dependencies_gradient_term_RHS(0, "grad(T)");
  set_dependencies_gradient_term_LHS(0, "grad(change(T))");

  set_variable_name(1, "q");
  set_variable_type(1, SCALAR);
  set_variable_equation_type(1, CONSTANT);

  set_variable_name(2, "error");
  set_variable_type(2, SCALAR);
  set_variable_equation_type(2, EXPLICIT_TIME_DEPENDENT);
  set_is_postprocessed_field(2, true);
  set_dependencies_value_term_RHS(2, "T");
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
      scalarGrad  Tx = variable_list.get_scalar_gradient(0);
      scalarValue q  = variable_list.get_scalar_value(1);

      variable_list.set_scalar_value_term(0, q);
      variable_list.set_scalar_gradient_term(0, -Tx);
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
      scalarGrad change_Tx = variable_list.get_scalar_gradient(0, CHANGE);

      variable_list.set_scalar_gradient_term(0, change_Tx, CHANGE);
    }
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_postprocess_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number> &variable_list,
  [[maybe_unused]] const dealii::Point<dim, dealii::VectorizedArray<number>> &q_point_loc)
  const
{
  scalarValue T = variable_list.get_scalar_value(0);

  scalarValue analytic =
    std::sin(M_PI * q_point_loc[0] / user_inputs.spatial_discretization.domain_size[0]) *
    q_point_loc[1] / user_inputs.spatial_discretization.domain_size[1] *
    (1.0 - q_point_loc[1] / user_inputs.spatial_discretization.domain_size[1]);

  scalarValue error = (T - analytic) * (T - analytic);

  variable_list.set_scalar_value_term(2, error);
}

INSTANTIATE_TRI_TEMPLATE(customPDE)

PRISMS_PF_END_NAMESPACE