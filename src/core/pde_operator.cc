// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/pde_operator.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
PDEOperator<dim, degree, number>::PDEOperator(
  const UserInputParameters<dim> &_user_inputs)
  : user_inputs(&_user_inputs)
{}

template <unsigned int dim, unsigned int degree, typename number>
const UserInputParameters<dim> &
PDEOperator<dim, degree, number>::get_user_inputs() const
{
  Assert(user_inputs != nullptr, dealii::ExcNotInitialized());
  return *user_inputs;
}

template <unsigned int dim, unsigned int degree, typename number>
number
PDEOperator<dim, degree, number>::get_timestep() const
{
  Assert(user_inputs != nullptr, dealii::ExcNotInitialized());
  return user_inputs->get_temporal_discretization().get_timestep();
}

#include "core/pde_operator.inst"

PRISMS_PF_END_NAMESPACE
