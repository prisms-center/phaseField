#ifndef UTILITIES_H
#define UTILITIES_H

#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

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
  for (unsigned int i = 0; i < input_vector.size(); ++i)
    {
      tensor[i] = dealii::make_vectorized_array(input_vector[i]);
    }

  return tensor;
}

#endif