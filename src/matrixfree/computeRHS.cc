// computeRHS() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"
#include "../../include/variableContainer.h"

// update RHS of each field
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::computeExplicitRHS()
{
  // log time
  computing_timer.enter_subsection("matrixFreePDE: computeRHS");

  // call to integrate and assemble while clearing residual vecotrs
  matrixFreeObject.cell_loop(&MatrixFreePDE<dim, degree>::getExplicitRHS,
                             this,
                             residualSet,
                             solutionSet,
                             true);

  // end log
  computing_timer.leave_subsection("matrixFreePDE: computeRHS");
}

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::getExplicitRHS(
  const MatrixFree<dim, double>               &data,
  std::vector<vectorType *>                   &dst,
  const std::vector<vectorType *>             &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  variableContainer<dim, degree, dealii::VectorizedArray<double>> variable_list(
    data,
    userInputs.varInfoListExplicitRHS);

  // loop over cells
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      variable_list.reinit_and_eval(src, cell);

      unsigned int num_q_points = variable_list.get_num_q_points();

      dealii::VectorizedArray<double> local_element_volume = element_volume[cell];

      // loop over quadrature points
      for (unsigned int q = 0; q < num_q_points; ++q)
        {
          variable_list.q_point = q;

          dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc =
            variable_list.get_q_point_location();

          // Calculate the residuals
          explicitEquationRHS(variable_list, q_point_loc, local_element_volume);
        }

      variable_list.integrate_and_distribute(dst);
    }
}

// update RHS of each field
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::computeNonexplicitRHS()
{
  // log time
  computing_timer.enter_subsection("matrixFreePDE: computeRHS");

  // call to integrate and assemble while clearing residual vecotrs
  matrixFreeObject.cell_loop(&MatrixFreePDE<dim, degree>::getNonexplicitRHS,
                             this,
                             residualSet,
                             solutionSet,
                             true);

  // end log
  computing_timer.leave_subsection("matrixFreePDE: computeRHS");
}

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::getNonexplicitRHS(
  const MatrixFree<dim, double>               &data,
  std::vector<vectorType *>                   &dst,
  const std::vector<vectorType *>             &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  variableContainer<dim, degree, dealii::VectorizedArray<double>> variable_list(
    data,
    userInputs.varInfoListNonexplicitRHS);

  // loop over cells
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      variable_list.reinit_and_eval(src, cell);

      unsigned int num_q_points = variable_list.get_num_q_points();

      dealii::VectorizedArray<double> local_element_volume = element_volume[cell];

      // loop over quadrature points
      for (unsigned int q = 0; q < num_q_points; ++q)
        {
          variable_list.q_point = q;

          dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc =
            variable_list.get_q_point_location();

          // Calculate the residuals
          nonExplicitEquationRHS(variable_list, q_point_loc, local_element_volume);
        }

      variable_list.integrate_and_distribute(dst);
    }
}

#include "../../include/matrixFreePDE_template_instantiations.h"
