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
customInitialCondition<dim>::set_initial_condition(
  [[maybe_unused]] const unsigned int             &index,
  [[maybe_unused]] const unsigned int             &component,
  [[maybe_unused]] const dealii::Point<dim>       &point,
  [[maybe_unused]] double                         &scalar_value,
  [[maybe_unused]] double                         &vector_component_value,
  [[maybe_unused]] const userInputParameters<dim> &user_inputs) const
{}

//JM helper function: a Kronecker Delta
int kDelta(unsigned int i, unsigned int j) 
{
  if (i == j) {
      return 1;
  } else {
      return 0;
  }
}

template <unsigned int dim, typename number>
void
customNonuniformDirichlet<dim, number>::set_nonuniform_dirichlet(
  [[maybe_unused]] const unsigned int               &index,
  [[maybe_unused]] const unsigned int               &boundary_id,
  [[maybe_unused]] const unsigned int               &component,
  [[maybe_unused]] const dealii::Point<dim>         &point,
  [[maybe_unused]] number                           &scalar_value,
  [[maybe_unused]] number                           &vector_component_value,
  [[maybe_unused]] const userInputParameters<dim>   &user_inputs) const
{
  
  const double radius = user_inputs.user_constants.get_model_constant_double("incRadius");
  const double poisson = user_inputs.user_constants.get_model_constant_double("poisson");
  const dealii::Tensor<1, dim, double> center =
  user_inputs.user_constants.get_model_constant_rank_1_tensor("center");

  double eshelbyConstant = (radius*radius*radius)/(6.0*(1-poisson));
  double dist_from_inclusion = 0.0;
  for (unsigned int i = 0; i < dim; i++)
    {
      //potential dist rewrite
      dist_from_inclusion = std::sqrt((point[0] - center[0]) * (point[0] - center[0]) +
                                      (point[1] - center[1]) * (point[1] - center[1]) +
                                      (point[2] - center[2]) * (point[2] - center[2]));

      double G = 0.0;
      for (unsigned int j = 0; j < dim; j ++)
        {
          for (unsigned int k = 0; k < dim; k++)
            {
              G += (kDelta(j,k)*0.01) * ((1-2*poisson) *
                   (kDelta(i,j)*((point[k]-center[k])/dist_from_inclusion) + 
                    kDelta(i,k)*((point[j]-center[j])/dist_from_inclusion) - 
                    kDelta(j,k)*((point[i]-center[i])/dist_from_inclusion)) + 
                    3 * ((point[i]-center[i])/dist_from_inclusion) * 
                    ((point[j]-center[j])/dist_from_inclusion) * 
                    ((point[k]-center[k])/dist_from_inclusion));
            }
        }
      
      if (component == i)
        {
          vector_component_value = -1.0*eshelbyConstant*(1/(dist_from_inclusion*dist_from_inclusion))*G;
        }
    }
}

template class customInitialCondition<1>;
template class customInitialCondition<2>;
template class customInitialCondition<3>;

template class customNonuniformDirichlet<1, double>;
template class customNonuniformDirichlet<2, double>;
template class customNonuniformDirichlet<3, double>;
template class customNonuniformDirichlet<1, float>;
template class customNonuniformDirichlet<2, float>;
template class customNonuniformDirichlet<3, float>;

PRISMS_PF_END_NAMESPACE
