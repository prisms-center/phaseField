#include <core/matrixFreePDE.h>

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::vmult(
  LinearAlgebra::distributed::Vector<double>       &dst,
  const LinearAlgebra::distributed::Vector<double> &src) const
{
  computing_timer.enter_subsection("matrixFreePDE: computeLHS");

  // Determine which function to use for the cell loop
  auto LHS_function = generatingInitialGuess ? &MatrixFreePDE<dim, degree>::getLaplaceLHS
                                             : &MatrixFreePDE<dim, degree>::getLHS;

  // Cell loop
  matrixFreeObject.cell_loop(LHS_function, this, dst, src, true);

  // Account for Dirichlet BC's (essentially copy dirichlet DOF values present in src to
  // dst, although it is unclear why the constraints can't just be distributed here)
  for (const auto &[index, value] : *valuesDirichletSet[currentFieldIndex])
    {
      if (dst.in_local_range(index))
        {
          dst(index) = src(index);
        }
    }

  computing_timer.leave_subsection("matrixFreePDE: computeLHS");
}

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::getLHS(
  const MatrixFree<dim, double>                    &data,
  LinearAlgebra::distributed::Vector<double>       &dst,
  const LinearAlgebra::distributed::Vector<double> &src,
  const std::pair<unsigned int, unsigned int>      &cell_range) const
{
  variableContainer<dim, degree, double> variable_list(
    data,
    userInputs.variable_info_list_LHS,
    userInputs.old_variable_info_list_LHS,
    userInputs.change_variable_info_list_LHS);

  // loop over cells
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      variable_list.reinit_and_eval(solutionSet, cell);
      variable_list.reinit_and_eval_change_in_solution(src, cell, currentFieldIndex);

      unsigned int num_q_points = variable_list.get_num_q_points();

      dealii::VectorizedArray<double> local_element_volume = element_volume[cell];

      // loop over quadrature points
      for (unsigned int q = 0; q < num_q_points; ++q)
        {
          variable_list.q_point = q;

          dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc =
            variable_list.get_q_point_location();

          // Calculate the residuals
          equationLHS(variable_list, q_point_loc, local_element_volume);
        }

      // Integrate the residuals and distribute from local to global
      variable_list.integrate_and_distribute_change_in_solution_LHS(dst,
                                                                    currentFieldIndex);
    }
}
