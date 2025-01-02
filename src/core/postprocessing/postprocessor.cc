#include <core/matrixFreePDE.h>

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::computePostProcessedFields(
  std::vector<vectorType *> &postProcessedSet)
{
  // Initialize the postProcessedSet
  for (unsigned int fieldIndex = 0; fieldIndex < pp_attributes.size(); fieldIndex++)
    {
      vectorType *U = new vectorType;
      postProcessedSet.push_back(U);
      matrixFreeObject.initialize_dof_vector(*U, 0);
    }

  // call to integrate and assemble
  matrixFreeObject.cell_loop(&MatrixFreePDE<dim, degree>::getPostProcessedFields,
                             this,
                             postProcessedSet,
                             solutionSet,
                             true);
}

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::getPostProcessedFields(
  const dealii::MatrixFree<dim, double>       &data,
  std::vector<vectorType *>                   &dst,
  const std::vector<vectorType *>             &src,
  const std::pair<unsigned int, unsigned int> &cell_range)
{
  // initialize FEEvaulation objects
  variableContainer<dim, degree, dealii::VectorizedArray<double>>
    variable_list(data, var_attributes, solveType::POSTPROCESS);
  variableContainer<dim, degree, dealii::VectorizedArray<double>>
    pp_variable_list(data, pp_attributes, solveType::POSTPROCESS);
  // loop over cells
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Initialize, read DOFs, and set evaulation flags for each variable
      variable_list.reinit_and_eval(src, cell);
      pp_variable_list.reinit(cell);

      unsigned int num_q_points = variable_list.get_num_q_points();

      dealii::VectorizedArray<double> local_element_volume = element_volume[cell];

      // loop over quadrature points
      for (unsigned int q = 0; q < num_q_points; ++q)
        {
          variable_list.q_point    = q;
          pp_variable_list.q_point = q;

          dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc =
            variable_list.get_q_point_location();

          // Calculate the residuals
          postProcessedFields(variable_list,
                              pp_variable_list,
                              q_point_loc,
                              local_element_volume);
        }

      pp_variable_list.integrate_and_distribute(dst);
    }
}
