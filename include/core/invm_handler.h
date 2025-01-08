#ifndef invm_handler_h
#define invm_handler_h

#include <deal.II/lac/la_parallel_vector.h>

#include <core/exceptions.h>
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
  ~invmHandler();

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
   * \brief Variable attributes. This is used to determine the proper return type for the
   * invm when given a field index.
   */
  const AttributesList &variable_attributes;

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
};

template <int dim, int degree>
invmHandler<dim, degree>::invmHandler(const AttributesList &_variable_attributes)
  : variable_attributes(_variable_attributes)
{
  for (const auto &[index, variable] : variable_attributes)
    {
      if (variable.field_type == fieldType::SCALAR)
        {
          scalar_needed = true;
        }
      if (variable.field_type == fieldType::VECTOR)
        {
          vector_needed = true;
        }
    }
}

template <int dim, int degree>
inline void
invmHandler<dim, degree>::compute_invm()
{}

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
  Assert(variable_attributes.at(index) != variable_attributes.end(),
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

  Assert(false, dealii::ExcMessage("Invalid field type."));
}

#endif