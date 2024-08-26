/*
 * vectorBCFunction.h
 *
 *  Created on: Feb 22, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_VECTORBCFUNCTION_H_
#define INCLUDE_VECTORBCFUNCTION_H_

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

template <int dim>
class vectorBCFunction : public dealii::Function<dim, double>
{
public:
  vectorBCFunction(const std::vector<double> BC_values);
  virtual void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<double> &values) const;

  virtual void
  vector_value_list(const std::vector<dealii::Point<dim>> &points,
                    std::vector<dealii::Vector<double>>   &value_list) const;

private:
  const std::vector<double> BC_values;
};

#endif /* INCLUDE_VECTORBCFUNCTION_H_ */
