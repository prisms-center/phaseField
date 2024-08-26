#ifndef INCLUDE_NONUNIFORMDIRICHLETBCS_H_
#define INCLUDE_NONUNIFORMDIRICHLETBCS_H_

#include "matrixFreePDE.h"
#include "userInputParameters.h"

template <int dim, int degree>
class NonUniformDirichletBC : public dealii::Function<dim>
{
public:
  dealii::Vector<double> values;

  NonUniformDirichletBC(const unsigned int          _index,
                        const unsigned int          _direction,
                        const double                _time,
                        MatrixFreePDE<dim, degree> *_matrix_free_pde)
    : dealii::Function<dim>(1)
    , index(_index)
    , direction(_direction)
    , time(_time)
    , matrix_free_pde(_matrix_free_pde)
  {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) + 1);
  }

  // IC for scalar values
  double
  value(const dealii::Point<dim> &p, const unsigned int component = 0) const
  {
    double                 scalar_BC = 0.0;
    dealii::Vector<double> vector_BC(dim);

    matrix_free_pde
      ->setNonUniformDirichletBCs(p, index, direction, time, scalar_BC, vector_BC);

    return scalar_BC;
  };

private:
  const unsigned int          index;
  const unsigned int          direction;
  const double                time;
  MatrixFreePDE<dim, degree> *matrix_free_pde;
};

template <int dim, int degree>
class NonUniformDirichletBCVector : public dealii::Function<dim>
{
public:
  dealii::Vector<double> values;

  NonUniformDirichletBCVector(const unsigned int          _index,
                              const unsigned int          _direction,
                              const double                _time,
                              MatrixFreePDE<dim, degree> *_matrix_free_pde)
    : dealii::Function<dim>(dim)
    , index(_index)
    , direction(_direction)
    , time(_time)
    , matrix_free_pde(_matrix_free_pde)
  {
    std::srand(dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) + 1);
  }

  // IC for vector values
  void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<double> &vector_BC) const
  {
    double scalar_BC = 0.0;

    vector_BC.reinit(dim);
    matrix_free_pde
      ->setNonUniformDirichletBCs(p, index, direction, time, scalar_BC, vector_BC);
  };

private:
  const unsigned int          index;
  const unsigned int          direction;
  const double                time;
  MatrixFreePDE<dim, degree> *matrix_free_pde;
};

#endif // INCLUDE_NONUNIFORMDIRICHLETBCS_H_
