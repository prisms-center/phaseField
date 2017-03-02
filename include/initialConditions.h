/*
 * initialConditions.h
 *
 *  Created on: Feb 27, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_INITIALCONDITIONS_H_
#define INCLUDE_INITIALCONDITIONS_H_

template <int dim>
class InitialCondition : public dealii::Function<dim>
{
public:
  unsigned int index;
  dealii::Vector<double> values;
  InitialCondition (const unsigned int _index) : dealii::Function<dim>(1), index(_index) {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const dealii::Point<dim> &p, const unsigned int component=0) const;

};

template <int dim>
class InitialConditionVec : public dealii::Function<dim>
{
public:
  unsigned int index;
  //Vector<double> values;
  InitialConditionVec (const unsigned int _index) : dealii::Function<dim>(dim), index(_index) {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  void vector_value (const dealii::Point<dim> &p,dealii::Vector<double> &vector_IC) const;

};

#endif /* INCLUDE_INITIALCONDITIONS_H_ */
