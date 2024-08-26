
#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::computeIntegral(double                   &integratedField,
                                            int                       index,
                                            std::vector<vectorType *> variableSet)
{
  QGauss<dim>         quadrature_formula(degree + 1);
  FE_Q<dim>           FE(QGaussLobatto<1>(degree + 1));
  FEValues<dim>       fe_values(FE,
                          quadrature_formula,
                          update_values | update_JxW_values | update_quadrature_points);
  const unsigned int  n_q_points = quadrature_formula.size();
  std::vector<double> cVal(n_q_points);

  // constraintsDirichletSet[index]->distribute(*variableSet[index]);
  // constraintsOtherSet[index]->distribute(*variableSet[index]);
  // variableSet[index]->update_ghost_values();

  typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandlersSet[0]
                                                          ->begin_active(),
                                                 endc = this->dofHandlersSet[0]->end();

  double value = 0.0;

  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          fe_values.get_function_values(*variableSet[index], cVal);

          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              value += (cVal[q]) * fe_values.JxW(q);
            }
        }
    }

  value = Utilities::MPI::sum(value, MPI_COMM_WORLD);

  // if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0){
  // std::cout<<"Integrated field: "<<value<<std::endl;
  // }

  integratedField = value;
}

//-----------------------------------

// update RHS of each field
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::computeIntegralMF(double &integratedField,
                                              int     index,
                                              const std::vector<vectorType *> variableSet)
{
  // log time
  computing_timer.enter_subsection("matrixFreePDE: computeIntegralMF");

  integrated_var = 0.0;
  integral_index = index;

  // call to integrate and assemble
  matrixFreeObject.cell_loop(&MatrixFreePDE<dim, degree>::getIntegralMF,
                             this,
                             residualSet,
                             variableSet);

  integratedField = Utilities::MPI::sum(integrated_var, MPI_COMM_WORLD);

  // end log
  computing_timer.leave_subsection("matrixFreePDE: computeIntegralMF");
}

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::getIntegralMF(
  const MatrixFree<dim, double>               &data,
  std::vector<vectorType *>                   &dst,
  const std::vector<vectorType *>             &src,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  dealii::FEEvaluation<dim, degree, degree + 1, 1, double> var(data, 0);

  // loop over cells
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      var.reinit(cell);
      var.read_dof_values_plain(*src[0]);
      var.evaluate(true, false, false);

      unsigned int num_q_points = var.n_q_points;

      // loop over quadrature points
      for (unsigned int q = 0; q < num_q_points; ++q)
        {
          dealii::VectorizedArray<double> val = var.get_value(q);
          dealii::VectorizedArray<double> jxw = var.JxW(q);
          assembler_lock.lock();
          for (unsigned i = 0; i < val.size(); i++)
            {
              integrated_var += val[i] * jxw[i];
            }
          assembler_lock.unlock();
        }
    }
}

#include "../../include/matrixFreePDE_template_instantiations.h"
