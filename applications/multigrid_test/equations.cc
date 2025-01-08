#include "custom_pde.h"

#include <config.h>
#include <core/type_enums.h>
#include <core/variable_attribute_loader.h>

void
customAttributeLoader::loadVariableAttributes()
{
  set_variable_name(0, "phi");
  set_variable_type(0, SCALAR);
  set_variable_equation_type(0, TIME_INDEPENDENT);

  set_dependencies_value_term_RHS(0, "");
  set_dependencies_gradient_term_RHS(0, "grad(phi)");
  set_dependencies_value_term_LHS(0, "");
  set_dependencies_gradient_term_LHS(0, "grad(change(phi))");
}

void
customAttributeLoader::loadPostProcessorVariableAttributes()
{}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_explicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number>    &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<number>> &q_point_loc) const
{}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_RHS(
  [[maybe_unused]] variableContainer<dim, degree, number>    &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<number>> &q_point_loc) const
{
  scalarGrad phix = variable_list.get_scalar_gradient(0);

  variable_list.set_scalar_gradient_term(0, -phix);
}

template <int dim, int degree, typename number>
void
customPDE<dim, degree, number>::compute_nonexplicit_LHS(
  [[maybe_unused]] variableContainer<dim, degree, number>    &variable_list,
  [[maybe_unused]] const Point<dim, VectorizedArray<number>> &q_point_loc) const
{
  scalarGrad change_phix = variable_list.get_scalar_gradient(0, CHANGE);

  variable_list.set_scalar_gradient_term(0, change_phix, CHANGE);
}

INSTANTIATE_TRI_TEMPLATE(customPDE)