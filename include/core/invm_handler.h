#ifndef invm_handler_h
#define invm_handler_h

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <core/exceptions.h>
#include <core/type_enums.h>
#include <core/user_inputs/user_input_parameters.h>
#include <core/variable_attributes.h>

/**
 * \brief This class handles the computation and access of the inverted mass matrix for
 * explicit solves.
 */
template <int dim, int degree>
class invmHandler
{
public:
  /**
   * \brief Constructor.
   */
  invmHandler(const AttributesList &_variable_attributes);

  /**
   * \brief Destructor.
   */
  ~invmHandler() = default;

  /**
   * \brief Initialize.
   */
  void
  initialize(std::shared_ptr<dealii::MatrixFree<dim, double>> _data);

  /**
   * \brief Compute the mass matrix for scalar/vector fields.
   */
  void
  compute_invm();

  /**
   * \brief Recompute the mass matrix for scalar/vector fields. This just points to
   * compute_invm() and is used for style.
   */
  void
  recompute_invm();

  /**
   * \brief  Getter function for the mass matrix for the given field index (constant
   * reference).
   */
  const dealii::LinearAlgebra::distributed::Vector<double> &
  get_invm(const uint &index);

private:
  /**
   * \brief Local computation of inverted mass matrix.
   */
  void
  compute_local_invm(const dealii::MatrixFree<dim, double>              &data,
                     dealii::LinearAlgebra::distributed::Vector<double> &dst,
                     const fieldType                                    &field_type,
                     const std::pair<uint, uint> &cell_range) const;

  /**
   * \brief Variable attributes. This is used to determine the proper return type for the
   * invm when given a field index.
   */
  const AttributesList &variable_attributes;

  /**
   * \brief Matrix-free object.
   */
  std::shared_ptr<dealii::MatrixFree<dim, double>> data;

  /**
   * \brief Inverse of the mass matrix for scalar fields.
   */
  dealii::LinearAlgebra::distributed::Vector<double> invm_scalar;

  /**
   * \brief Inverse of the mass matrix for vector fields.
   */
  dealii::LinearAlgebra::distributed::Vector<double> invm_vector;

  /**
   * \brief Whether a scalar invm is needed.
   */
  bool scalar_needed = false;

  /**
   * \brief Whether a vector invm is needed.
   */
  bool vector_needed = false;

  /**
   * \brief Field index of the first occuring scalar field. This is the index for which we
   * attached the FEEvaluation objects to evaluate and initialize the invm vector.
   */
  uint scalar_index;

  /**
   * \brief Field index of the first occuring vector field. This is the index for which we
   * attached the FEEvaluation objects to evaluate and initialize the invm vector.
   */
  uint vector_index;
};

template <int dim, int degree>
invmHandler<dim, degree>::invmHandler(const AttributesList &_variable_attributes)
  : variable_attributes(_variable_attributes)
{
  for (const auto &[index, variable] : variable_attributes)
    {
      if (variable.field_type == fieldType::SCALAR && !scalar_needed)
        {
          scalar_needed = true;
          scalar_index  = index;
        }
      if (variable.field_type == fieldType::VECTOR && !vector_needed)
        {
          vector_needed = true;
          vector_index  = index;
        }
    }
}

template <int dim, int degree>
inline void
invmHandler<dim, degree>::initialize(
  std::shared_ptr<dealii::MatrixFree<dim, double>> _data)
{
  // Grab the shared_ptr to the matrix-free object
  data = _data;
}

template <int dim, int degree>
inline void
invmHandler<dim, degree>::compute_invm()
{
  // Initialize the invm vectors and cell loop to compute the invm vector, as neccessary
  if (scalar_needed)
    {
      data->initialize_dof_vector(invm_scalar, scalar_index);
      data->cell_loop(&invmHandler::compute_local_invm,
                      this,
                      invm_scalar,
                      fieldType::SCALAR);

      // Loop over cells and take the inverse
      for (uint i = 0; i < invm_scalar.locally_owned_size(); ++i)
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
      data->cell_loop(&invmHandler::compute_local_invm,
                      this,
                      invm_vector,
                      fieldType::VECTOR);

      // Loop over cells and take the inverse
      for (uint i = 0; i < invm_vector.locally_owned_size(); ++i)
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
inline void
invmHandler<dim, degree>::recompute_invm()
{
  this->compute_invm();
}

template <int dim, int degree>
inline const dealii::LinearAlgebra::distributed::Vector<double> &
invmHandler<dim, degree>::get_invm(const uint &index)
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
               "variable attributes correspond with the provided index."));
      Assert(invm_scalar.size() == 0,
             dealii::ExcMessage("The scalar invm has size 0. Please make sure to call "
                                "compute_invm() prior to calling the getter function."));

      return invm_scalar;
    }
  else if (variable_attributes.at(index).field_type == fieldType::VECTOR)
    {
      Assert(vector_needed,
             dealii::ExcMessage(
               "The invm for vector fields is marked as not needed. Make sure the "
               "variable attributes correspond with the provided index."));
      Assert(invm_vector.size() == 0,
             dealii::ExcMessage("The vector invm has size 0. Please make sure to call "
                                "compute_invm() prior to calling the getter function."));

      return invm_vector;
    }

  Assert(false, dealii::ExcMessage("Invalid field type"));
  throw std::runtime_error("Invalid field type");
}

template <int dim, int degree>
inline void
invmHandler<dim, degree>::compute_local_invm(
  const dealii::MatrixFree<dim, double>              &data,
  dealii::LinearAlgebra::distributed::Vector<double> &dst,
  const fieldType                                    &field_type,
  const std::pair<uint, uint>                        &cell_range) const
{
  if (field_type == fieldType::SCALAR)
    {
      dealii::FEEvaluation<dim, degree, degree + 1, 1, double> fe_eval(data,
                                                                       scalar_index);

      for (uint cell = cell_range.first; cell < cell_range.second; ++cell)
        {
          fe_eval.reinit(cell);
          for (const uint q : fe_eval.quadrature_point_indices())
            {
              fe_eval.submit_value(dealii::VectorizedArray<double>(1.0), q);
            }
          fe_eval.integrate(dealii::EvaluationFlags::values);
          fe_eval.distribute_local_to_global(dst);
        }
    }
  else if (field_type == fieldType::VECTOR)
    {
      dealii::FEEvaluation<dim, degree, degree + 1, dim, double> fe_eval(data,
                                                                         vector_index);

      dealii::Tensor<1, dim, dealii::VectorizedArray<double>> one;
      one = 1.0;

      for (uint cell = cell_range.first; cell < cell_range.second; ++cell)
        {
          fe_eval.reinit(cell);
          for (const uint q : fe_eval.quadrature_point_indices())
            {
              fe_eval.submit_value(one, q);
            }
          fe_eval.integrate(dealii::EvaluationFlags::values);
          fe_eval.distribute_local_to_global(dst);
        }
    }
}

#endif