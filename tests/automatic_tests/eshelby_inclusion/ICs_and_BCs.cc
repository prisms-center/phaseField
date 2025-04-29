// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/nonuniform_dirichlet.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

#include <cmath>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
void
customInitialCondition<dim>::setInitialCondition(
  [[maybe_unused]] const dealii::Point<dim>       &point,
  [[maybe_unused]] const unsigned int             &index,
  [[maybe_unused]] const unsigned int             &component,
  [[maybe_unused]] double                         &scalar_value,
  [[maybe_unused]] double                         &vector_component_value,
  [[maybe_unused]] const userInputParameters<dim> &user_inputs) const
{}

//JM helper function: a Kronecker Delta
int kDelta(int i, int j) 
{
  if (i == j) {
      return 1;
  } else {
      return 0;
  }
}

template <int dim, int degree>
void
customNonuniformDirichlet<dim, degree>::set_nonuniform_dirichlet(
  [[maybe_unused]] const dealii::Point<dim>         &point,
  [[maybe_unused]] const unsigned int               &index,
  [[maybe_unused]] const unsigned int               &boundary_id,
  [[maybe_unused]] const unsigned int               &component,
  [[maybe_unused]] number                           &vector_component_value,
  [[maybe_unused]] number                           &scalar_value,
  [[maybe_unused]] const userInputerParameters<dim> &user_inputs)
{
  //JM need to update this BC
  for (unsigned int i = 0; i < dim; i++)
  {
    scalarValue eshelbyConstant = (incRadius*incradius*incRadius)/(6.0*(1-poisson));
    scalarValue dist_from_inclusion
    for (unsigned int i = 0; i < dim; i++)
      {
        dist_from_inclusion += (point[i] - constV<number>(0.0)) *
                               (point[i] - constV<number>(0.0));
      }
    scalarValue G;
    for (unsigned int j = 0; j < dim; j ++)
    {
      for (unsigned int k = 0; k < dim; k++)
      {
        G += (kDelta(j,k)*0.01) * ((1-2*poisson) *
             (kDelta(i,j)*((point[k]-center[k])/dist_from_inclusion) + 
              kDelta(i,k)*((point[j]-center[j])/dist_from_inclusion) - 
              kDelta(j,k)*((point[i]-center[i])/dist_from_inclusion)) + 
              3 * ((point[i]-center[k])/dist_from_inclusion) * 
              ((point[j]-center[j])/dist_from_inclusion) * 
              ((point[k]-center[i])/dist_from_inclusion));
      }
    }
    vector_BC(i) = -1.0*A*(1/(dist*dist))*G;
  }
}
