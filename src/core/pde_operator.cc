#include <prismspf/core/pde_operator.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, int degree, typename number>
PDEOperator<dim, degree, number>::PDEOperator(
  const userInputParameters<dim> &_user_inputs)
  : user_inputs(&_user_inputs)
{}

template <int dim, int degree, typename number>
const userInputParameters<dim> &
PDEOperator<dim, degree, number>::get_user_inputs() const
{
  Assert(user_inputs != nullptr, dealii::ExcNotInitialized());
  return *user_inputs;
}

template <int dim, int degree, typename number>
number
PDEOperator<dim, degree, number>::get_timestep() const
{
  Assert(user_inputs != nullptr, dealii::ExcNotInitialized());
  return user_inputs->temporal_discretization.dt;
}

INSTANTIATE_TRI_TEMPLATE(PDEOperator)

PRISMS_PF_END_NAMESPACE