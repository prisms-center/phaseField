/*
 * initialConditions.h
 *
 *  Created on: Feb 27, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_INITIALCONDITIONS_H_
#define INCLUDE_INITIALCONDITIONS_H_

#include "userInputParameters.h"
#include "matrixFreePDE.h"

template <int dim, int degree>
class InitialConditionScalar : public dealii::Function<dim>
{
public:
  const unsigned int index;
  const userInputParameters<dim> userInputs;
  dealii::Vector<double> values;
  InitialConditionScalar (const unsigned int _index, const userInputParameters<dim> _userInputs, MatrixFreePDE<dim,degree>* _mfPDE) : dealii::Function<dim>(1), index(_index), userInputs(_userInputs),mfPDE(_mfPDE) {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const dealii::Point<dim> &p, const unsigned int component=0) const {
      double scalar_IC = 0.0;
      dealii::Vector<double> vector_IC(dim);

      mfPDE->setInitialCondition(p, index, scalar_IC, vector_IC);
      return scalar_IC;
  };

  void vector_value (const dealii::Point<dim> &p,dealii::Vector<double> &vector_IC) const {
      double scalar_IC = 0.0;
      mfPDE->setInitialCondition(p, index, scalar_IC, vector_IC);

  };

private:
    MatrixFreePDE<dim,degree>* mfPDE;
};

/*
template <int dim, int degree>
class InitialConditionVector : public dealii::Function<dim>
{
public:
  const unsigned int index;
  const userInputParameters<dim> userInputs;
  dealii::Vector<double> values;
  InitialConditionVector (const unsigned int _index, const userInputParameters<dim> _userInputs, MatrixFreePDE<dim,degree>* _mfPDE) : dealii::Function<dim>(1), index(_index), userInputs(_userInputs),mfPDE(_mfPDE) {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }

  void vector_value (const dealii::Point<dim> &p,dealii::Vector<double> &vector_IC) const {
      double scalar_IC = 0.0;
      mfPDE->setInitialCondition(p, index, scalar_IC, vector_IC);

  };

private:
    MatrixFreePDE<dim,degree>* mfPDE;
};
*/

#endif /* INCLUDE_INITIALCONDITIONS_H_ */
