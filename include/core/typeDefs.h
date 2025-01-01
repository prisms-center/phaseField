#ifndef INCLUDE_TYPEDEFS_H_
#define INCLUDE_TYPEDEFS_H_

#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

template <int dim, typename number>
struct TypeDefs
{
  using scalarValue = dealii::VectorizedArray<number>;
  using scalarGrad  = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using scalarHess  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;

  using vectorValue = dealii::Tensor<1, dim, dealii::VectorizedArray<number>>;
  using vectorGrad  = dealii::Tensor<2, dim, dealii::VectorizedArray<number>>;
  using vectorHess  = dealii::Tensor<3, dim, dealii::VectorizedArray<number>>;
};

#endif
