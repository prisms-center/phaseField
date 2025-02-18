#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/config.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/user_inputs/user_input_parameters.h>

#include <map>
#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, int degree>
invmHandler<dim, degree>::invmHandler(
  const std::map<unsigned int, variableAttributes> &_variable_attributes)
  : variable_attributes(_variable_attributes)
{
  for (const auto &[index, variable] : variable_attributes)
    {
      if (variable.field_type == fieldType::SCALAR && !scalar_needed &&
          variable.pde_type == PDEType::EXPLICIT_TIME_DEPENDENT)
        {
          scalar_needed = true;
          scalar_index  = index;
        }
      if (variable.field_type == fieldType::VECTOR && !vector_needed &&
          variable.pde_type == PDEType::EXPLICIT_TIME_DEPENDENT)
        {
          vector_needed = true;
          vector_index  = index;
        }
    }
}

template <int dim, int degree>
void
invmHandler<dim, degree>::initialize(
  std::shared_ptr<dealii::MatrixFree<dim, double, size_type>> _data)
{
  // Grab the shared_ptr to the matrix-free object
  data = _data;
}

template <int dim, int degree>
void
invmHandler<dim, degree>::compute_invm()
{
  // Initialize the invm vectors and cell loop to compute the invm vector, as neccessary
  if (scalar_needed)
    {
      data->initialize_dof_vector(invm_scalar, scalar_index);

      dealii::FEEvaluation<dim, degree, degree + 1, 1, double> fe_eval(*data,
                                                                       scalar_index);

      for (unsigned int cell = 0; cell < data->n_cell_batches(); ++cell)
        {
          fe_eval.reinit(cell);
          for (const unsigned int q : fe_eval.quadrature_point_indices())
            {
              fe_eval.submit_value(dealii::VectorizedArray<double>(1.0), q);
            }
          fe_eval.integrate(dealii::EvaluationFlags::values);
          fe_eval.distribute_local_to_global(invm_scalar);
        }
      invm_scalar.compress(dealii::VectorOperation::add);

      // Loop over cells and take the inverse
      for (unsigned int i = 0; i < invm_scalar.locally_owned_size(); ++i)
        {
          if (invm_scalar.local_element(i) > 1.0e-15)
            {
              invm_scalar.local_element(i) = 1.0 / invm_scalar.local_element(i);
            }
          else
            {
              invm_scalar.local_element(i) = 1;
            }
        }
    }
  if (vector_needed)
    {
      data->initialize_dof_vector(invm_vector, vector_index);

      dealii::FEEvaluation<dim, degree, degree + 1, dim, double> fe_eval(*data,
                                                                         vector_index);

      dealii::Tensor<1, dim, dealii::VectorizedArray<double>> one;
      for (unsigned int i = 0; i < dim; i++)
        {
          one[i] = 1.0;
        }

      for (unsigned int cell = 0; cell < data->n_cell_batches(); ++cell)
        {
          fe_eval.reinit(cell);
          for (const unsigned int q : fe_eval.quadrature_point_indices())
            {
              fe_eval.submit_value(one, q);
            }
          fe_eval.integrate(dealii::EvaluationFlags::values);
          fe_eval.distribute_local_to_global(invm_vector);
        }
      invm_vector.compress(dealii::VectorOperation::add);

      // Loop over cells and take the inverse
      for (unsigned int i = 0; i < invm_vector.locally_owned_size(); ++i)
        {
          if (invm_vector.local_element(i) > 1.0e-15)
            {
              invm_vector.local_element(i) = 1.0 / invm_vector.local_element(i);
            }
          else
            {
              invm_vector.local_element(i) = 1;
            }
        }
    }
}

template <int dim, int degree>
void
invmHandler<dim, degree>::recompute_invm()
{
  this->compute_invm();
}

template <int dim, int degree>
const typename invmHandler<dim, degree>::VectorType &
invmHandler<dim, degree>::get_invm(const unsigned int &index) const
{
  Assert(variable_attributes.find(index) != variable_attributes.end(),
         dealii::ExcMessage(
           "Invalid index. The provided index does not have an entry in the variable "
           "attributes that were provided to the constructor."));

  if (variable_attributes.at(index).field_type == fieldType::SCALAR)
    {
      Assert(scalar_needed,
             dealii::ExcMessage(
               "The invm for scalar fields is marked as not needed. Make sure the "
               "variable attributes correspond with the provided index. Additionally, "
               "the invm is only necessary for explicit fields."));
      Assert(invm_scalar.size() != 0,
             dealii::ExcMessage("The scalar invm has size 0. Please make sure to call "
                                "compute_invm() prior to calling the getter function."));

      return invm_scalar;
    }
  else if (variable_attributes.at(index).field_type == fieldType::VECTOR)
    {
      Assert(vector_needed,
             dealii::ExcMessage(
               "The invm for vector fields is marked as not needed. Make sure the "
               "variable attributes correspond with the provided index. Additionally, "
               "the invm is only necessary for explicit fields."));
      Assert(invm_vector.size() != 0,
             dealii::ExcMessage("The vector invm has size 0. Please make sure to call "
                                "compute_invm() prior to calling the getter function."));

      return invm_vector;
    }

  AssertThrow(false, dealii::ExcMessage("Invalid field type"));
}

INSTANTIATE_BI_TEMPLATE(invmHandler)

PRISMS_PF_END_NAMESPACE