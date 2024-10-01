// computeInvM() method for MatrixFreePDE class
#include <deal.II/matrix_free/evaluation_flags.h>

#include "../../include/matrixFreePDE.h"
#include <functional>
#include <numeric>

// compute inverse of the diagonal mass matrix and store in vector invM
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::computeInvM()
{
  // Find invM indices by looping through all fields for the first field of that
  // type
  //(e.g., the index of the first explicit time dependent scalar field)
  // This has to be done for all types of fields.

  // This only has to be done when the fields are first initialized. It's
  // redundant for any reinitialized like in adaptive mesh refinement (AMR).
  // Hopefully the performance is fine
  bool invMscalarFound = false;
  bool invMvectorFound = false;

  unsigned int parabolicScalarFieldIndex = 0;
  unsigned int parabolicVectorFieldIndex = 0;

  for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
    {
      if (fields[fieldIndex].pdetype == EXPLICIT_TIME_DEPENDENT ||
          fields[fieldIndex].pdetype == AUXILIARY)
        {
          if (fields[fieldIndex].type == SCALAR && !invMscalarFound)
            {
              parabolicScalarFieldIndex = fieldIndex;
              invMscalarFound           = true;
              continue;
            }
          else if (fields[fieldIndex].type == VECTOR && !invMvectorFound)
            {
              parabolicVectorFieldIndex = fieldIndex;
              invMvectorFound           = true;
              continue;
            }
        }
    }

  // Check if invM has been found. If not, print an "error" message
  if (!invMscalarFound)
    {
      pcout << "matrixFreePDE.h: no PARABOLIC scalar field... hence setting "
               "parabolicScalarFieldIndex to 0 and marching ahead withn invM "
               "computation\n";
    }
  else if (!invMvectorFound)
    {
      pcout << "matrixFreePDE.h: no PARABOLIC vector field... hence setting "
               "parabolicVectorFieldIndex to 0 and marching ahead withn invM "
               "computation\n";
    }

  // Initialize invM and clear its values
  matrixFreeObject.initialize_dof_vector(invMscalar, parabolicScalarFieldIndex);
  invMscalar = 0.0;
  matrixFreeObject.initialize_dof_vector(invMvector, parabolicVectorFieldIndex);
  invMvector = 0.0;

  // invM evaluation flags
  dealii::EvaluationFlags::EvaluationFlags invM_flags = dealii::EvaluationFlags::values;

  // Compute mass matrix for the given type of quadrature. Selecting gauss
  // lobatto quadrature points which are suboptimal but give diagonal M
  if (fields[parabolicScalarFieldIndex].type == SCALAR)
    {
      VectorizedArray<double>   one = make_vectorized_array(1.0);
      FEEvaluation<dim, degree> fe_eval(matrixFreeObject, parabolicScalarFieldIndex);

      const unsigned int n_q_points = fe_eval.n_q_points;
      for (unsigned int cell = 0; cell < matrixFreeObject.n_cell_batches(); ++cell)
        {
          fe_eval.reinit(cell);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              fe_eval.submit_value(one, q);
            }
          fe_eval.integrate(invM_flags);
          fe_eval.distribute_local_to_global(invMscalar);
        }
    }
  if (fields[parabolicVectorFieldIndex].type == VECTOR)
    {
      dealii::Tensor<1, dim, dealii::VectorizedArray<double>> oneV;
      for (unsigned int i = 0; i < dim; i++)
        {
          oneV[i] = 1.0;
        }

      FEEvaluation<dim, degree, degree + 1, dim> fe_eval(matrixFreeObject,
                                                         parabolicVectorFieldIndex);

      const unsigned int n_q_points = fe_eval.n_q_points;
      for (unsigned int cell = 0; cell < matrixFreeObject.n_cell_batches(); ++cell)
        {
          fe_eval.reinit(cell);
          for (unsigned int q = 0; q < n_q_points; ++q)
            {
              fe_eval.submit_value(oneV, q);
            }
          fe_eval.integrate(invM_flags);
          fe_eval.distribute_local_to_global(invMvector);
        }
    }

  invMscalar.compress(VectorOperation::add);
  invMvector.compress(VectorOperation::add);

  // Calculate the volume of the smallest cell to prevent a non-zero value of
  // invM being confused for a near zero value (which can happen if the domain
  // size is 1e-6 or below)
  std::vector<double> min_element_length;
  for (unsigned int d = 0; d < dim; d++)
    {
      int num_elements =
        userInputs.subdivisions.at(d) *
        dealii::Utilities::fixed_power<2>(userInputs.max_refinement_level);
      min_element_length.push_back(userInputs.domain_size[d] / double(num_elements));
    }
  double min_cell_volume = std::accumulate(begin(min_element_length),
                                           end(min_element_length),
                                           1.0,
                                           std::multiplies<>());

  // Invert scalar mass matrix diagonal elements
#if (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR < 4)
  for (unsigned int k = 0; k < invMscalar.local_size(); ++k)
    {
#else
  for (unsigned int k = 0; k < invMscalar.locally_owned_size(); ++k)
    {
#endif
      if (std::abs(invMscalar.local_element(k)) > 1.0e-15 * min_cell_volume)
        {
          invMscalar.local_element(k) = 1. / invMscalar.local_element(k);
        }
      else
        {
          invMscalar.local_element(k) = 0;
        }
    }
  pcout << "computed scalar mass matrix (using FE space for field: "
        << parabolicScalarFieldIndex << ")\n";

  // Invert vector mass matrix diagonal elements
#if (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR < 4)
  for (unsigned int k = 0; k < invMvector.local_size(); ++k)
    {
#else
  for (unsigned int k = 0; k < invMvector.locally_owned_size(); ++k)
    {
#endif
      if (std::abs(invMvector.local_element(k)) > 1.0e-15 * min_cell_volume)
        {
          invMvector.local_element(k) = 1. / invMvector.local_element(k);
        }
      else
        {
          invMvector.local_element(k) = 0;
        }
    }
  pcout << "computed vector mass matrix (using FE space for field: "
        << parabolicVectorFieldIndex << ")\n";
}

#include "../../include/matrixFreePDE_template_instantiations.h"
