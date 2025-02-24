// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

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

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, int degree, typename number>
invmHandler<dim, degree, number>::invmHandler(
  const std::map<unsigned int, variableAttributes> &_variable_attributes)
  : variable_attributes(_variable_attributes)
  , data(nullptr)
  , invm_scalar(VectorType())
  , invm_vector(VectorType())
{
  for (const auto &[index, variable] : variable_attributes)
    {
      if (variable.field_type == fieldType::SCALAR && !scalar_needed &&
          (variable.pde_type == PDEType::EXPLICIT_TIME_DEPENDENT ||
           variable.pde_type == PDEType::AUXILIARY))
        {
          scalar_needed = true;
          scalar_index  = index;
        }
      if (variable.field_type == fieldType::VECTOR && !vector_needed &&
          (variable.pde_type == PDEType::EXPLICIT_TIME_DEPENDENT ||
           variable.pde_type == PDEType::AUXILIARY))
        {
          vector_needed = true;
          vector_index  = index;
        }
    }
}

template <int dim, int degree, typename number>
void
invmHandler<dim, degree, number>::initialize(
  std::shared_ptr<dealii::MatrixFree<dim, number, size_type>> _data)
{
  Assert(data == nullptr,
         dealii::ExcMessage("A ptr to a matrix-free object has already been assigned. "
                            "Please call clear() before calling this function."));

  // Grab the shared_ptr to the matrix-free object
  data = _data;
}

template <int dim, int degree, typename number>
void
invmHandler<dim, degree, number>::compute_invm()
{
  Assert(data != nullptr, dealii::ExcNotInitialized());

  // Initialize the invm vectors and cell loop to compute the invm vector, as neccessary
  if (scalar_needed)
    {
      data->initialize_dof_vector(invm_scalar, scalar_index);

      dealii::FEEvaluation<dim, degree, degree + 1, 1, number> fe_eval(*data,
                                                                       scalar_index);

      for (unsigned int cell = 0; cell < data->n_cell_batches(); ++cell)
        {
          fe_eval.reinit(cell);
          for (const unsigned int q : fe_eval.quadrature_point_indices())
            {
              fe_eval.submit_value(dealii::VectorizedArray<number>(1.0), q);
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

      dealii::FEEvaluation<dim, degree, degree + 1, dim, number> fe_eval(*data,
                                                                         vector_index);

      dealii::Tensor<1, dim, dealii::VectorizedArray<number>> one;
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

template <int dim, int degree, typename number>
void
invmHandler<dim, degree, number>::recompute_invm()
{
  this->compute_invm();
}

template <int dim, int degree, typename number>
const typename invmHandler<dim, degree, number>::VectorType &
invmHandler<dim, degree, number>::get_invm(const unsigned int &index) const
{
  Assert(data != nullptr, dealii::ExcNotInitialized());

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

  AssertThrow(false, UnreachableCode());
}

template <int dim, int degree, typename number>
void
invmHandler<dim, degree, number>::clear()
{
  data          = nullptr;
  invm_scalar   = VectorType();
  invm_vector   = VectorType();
  scalar_needed = false;
  vector_needed = false;
  scalar_index  = numbers::invalid_index;
  vector_index  = numbers::invalid_index;
}

INSTANTIATE_TRI_TEMPLATE(invmHandler)

PRISMS_PF_END_NAMESPACE