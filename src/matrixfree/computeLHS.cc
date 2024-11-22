// vmult() and getLHS() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

// vmult operation for LHS
template <int dim, int degree>

void
MatrixFreePDE<dim, degree>::vmult(vectorType &dst, const vectorType &src) const
{
  // log time
  computing_timer.enter_subsection("matrixFreePDE: computeLHS");

  // create temporary copy of src vector as src2, as vector src is marked const
  // and cannot be changed
  dealii::LinearAlgebra::distributed::Vector<double> src2;
  matrixFreeObject.initialize_dof_vector(src2, currentFieldIndex);
  src2 = src;

  // call cell_loop
  if (!generatingInitialGuess)
    {
      matrixFreeObject.cell_loop(&MatrixFreePDE<dim, degree>::getLHS,
                                 this,
                                 dst,
                                 src2,
                                 true);
    }
  else
    {
      matrixFreeObject.cell_loop(&MatrixFreePDE<dim, degree>::getLaplaceLHS,
                                 this,
                                 dst,
                                 src2,
                                 true);
    }

  // Account for Dirichlet BC's (essentially copy dirichlet DOF values present in src to
  // dst, although it is unclear why the constraints can't just be distributed here)
  for (auto &it : *valuesDirichletSet[currentFieldIndex])
    {
      if (dst.in_local_range(it.first))
        {
          dst(it.first) = src(it.first); //*jacobianDiagonal(it->first);
        }
    }

  // end log
  computing_timer.leave_subsection("matrixFreePDE: computeLHS");
}

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::getLHS(
  const MatrixFree<dim, double>               &data,
  vectorType                                  &dst,
  const vectorType                            &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  variableContainer<dim, degree, dealii::VectorizedArray<double>>
    variable_list(data, userInputs.varInfoListLHS, userInputs.varChangeInfoListLHS);

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

#include "../../include/matrixFreePDE_template_instantiations.h"
