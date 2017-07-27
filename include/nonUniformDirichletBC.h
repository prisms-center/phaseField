#ifndef INCLUDE_NONUNIFORMDIRICHLETBCS_H_
#define INCLUDE_NONUNIFORMDIRICHLETBCS_H_

#include "userInputParameters.h"

template <int dim>
class NonUniformDirichletBC : public dealii::Function<dim>
{
public:
  const unsigned int index;
  const unsigned int direction;

  const userInputParameters<dim> userInputs;

  NonUniformDirichletBC (const unsigned int _index, const unsigned int _direction, const userInputParameters<dim> _userInputs) : dealii::Function<dim>(1), index(_index), direction(_direction), userInputs(_userInputs) {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  double value (const dealii::Point<dim> &p, const unsigned int component=0) const;

 };

 template <int dim>
 class NonUniformDirichletBCVec : public dealii::Function<dim>
 {
 public:
   const unsigned int index;
   const unsigned int direction;

   const userInputParameters<dim> userInputs;

   NonUniformDirichletBCVec (const unsigned int _index, const unsigned int _direction, const userInputParameters<dim> _userInputs) : dealii::Function<dim>(1), index(_index), direction(_direction), userInputs(_userInputs) {
     std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
   }
   void vector_value (const dealii::Point<dim> &p,dealii::Vector<double> &vector_BC) const;

  };


#endif // INCLUDE_NONUNIFORMDIRICHLETBCS_H_
