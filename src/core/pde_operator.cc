// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/pde_operator.h>
#include <prismspf/core/phase_field_tools.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
PDEOperator<dim, degree, number>::PDEOperator(
  const UserInputParameters<dim> &_user_inputs,
  const PhaseFieldTools<dim>     &_pf_tools)
  : user_inputs(&_user_inputs)
  , pf_tools(&_pf_tools)
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
  return static_cast<number>(user_inputs->get_temporal_discretization().get_timestep());
}

template <unsigned int dim, unsigned int degree, typename number>
const PhaseFieldTools<dim> &
PDEOperator<dim, degree, number>::get_pf_tools() const
{
  Assert(pf_tools != nullptr, dealii::ExcNotInitialized());
  return *pf_tools;
}

#include "core/pde_operator.inst"

PRISMS_PF_END_NAMESPACE
