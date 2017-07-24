/*
 * initialConditions.h
 *
 *  Created on: Feb 27, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_INITIALCONDITIONS_H_
#define INCLUDE_INITIALCONDITIONS_H_

#include "userInputParameters.h"

template <int dim>
class InitialCondition : public dealii::Function<dim>
{
public:
  const unsigned int index;
  const userInputParameters<dim> userInputs;
  dealii::Vector<double> values;
  InitialCondition (const unsigned int _index, const userInputParameters<dim> _userInputs) : dealii::Function<dim>(1), index(_index), userInputs(_userInputs) {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const dealii::Point<dim> &p, const unsigned int component=0) const;

};

template <int dim>
class InitialConditionVec : public dealii::Function<dim>
{
public:
  const unsigned int index;
  const userInputParameters<dim> userInputs;
  //Vector<double> values;
  InitialConditionVec (const unsigned int _index, const userInputParameters<dim> _userInputs) : dealii::Function<dim>(dim), index(_index), userInputs(_userInputs) {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  void vector_value (const dealii::Point<dim> &p,dealii::Vector<double> &vector_IC) const;

};

#endif /* INCLUDE_INITIALCONDITIONS_H_ */
