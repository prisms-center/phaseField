/*
 * vectorBCFunction.cc
 *
 *  Created on: Feb 22, 2017
 *      Author: stephendewitt
 */

#include "../../include/vectorBCFunction.h"

#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>

template <int dim>
vectorBCFunction<dim>::vectorBCFunction(const std::vector<double> &BC_values)
  : dealii::Function<dim>(dim)
  , BC_values(BC_values)
{}

template <int dim>
void
vectorBCFunction<dim>::vector_value([[maybe_unused]] const dealii::Point<dim> &p,
                                    dealii::Vector<double> &values) const
{
  for (unsigned int i = 0; i < dim; i++)
    {
      values(i) = BC_values[i];
    }
}

template <int dim>
void
vectorBCFunction<dim>::vector_value_list(
  const std::vector<dealii::Point<dim>> &points,
  std::vector<dealii::Vector<double>>   &value_list) const
{
  const unsigned int n_points = points.size();
  for (unsigned int p = 0; p < n_points; ++p)
    {
      vectorBCFunction<dim>::vector_value(points[p], value_list[p]);
    }
}

template class vectorBCFunction<1>;
template class vectorBCFunction<2>;
template class vectorBCFunction<3>;
