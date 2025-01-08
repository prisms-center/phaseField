#ifndef utilities_h
#define utilities_h

#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/evaluation_flags.h>

#include <core/exceptions.h>
#include <vector>

/**
 * \brief Convert given scalar to vectorized array.
 */
#define constV(value) make_vectorized_array(value)

/**
 * \brief Convert given vector to a tensor of vectorized arrays.
 */
#define constT(vector, dim) make_tensor_of_vectorized_arrays<dim>(vector)

template <int dim, typename datatype>
dealii::Tensor<1, dim, dealii::VectorizedArray<datatype>>
make_tensor_of_vectorized_arrays(const std::vector<datatype> &input_vector)
{
  AssertDimension(input_vector.size(), dim);

  dealii::Tensor<1, dim, dealii::VectorizedArray<datatype>> tensor;

  // Populate the Tensor with vectorized arrays
  for (uint i = 0; i < input_vector.size(); ++i)
    {
      tensor[i] = dealii::make_vectorized_array(input_vector[i]);
    }

  return tensor;
}

/**
 * \brief Convert bool to string.
 */
inline const char *
bool_to_string(bool b)
{
  return b ? "true" : "false";
}

/**
 * \brief Convert evaluation flags to string.
 */
inline std::string
eval_flags_to_string(dealii::EvaluationFlags::EvaluationFlags flag)
{
  std::string result;

  if (flag & dealii::EvaluationFlags::values)
    {
      result += "values";
    }
  if (flag & dealii::EvaluationFlags::gradients)
    {
      if (!result.empty())
        {
          result += " | ";
        }
      result += "gradients";
    }
  if (flag & dealii::EvaluationFlags::hessians)
    {
      if (!result.empty())
        {
          result += " | ";
        }
      result += "hessians";
    }

  if (result.empty())
    {
      return "nothing";
    }

  return result;
}

#endif