// Model Variables Class

#ifndef INCLUDE_MODELVARIABLE_H_
#define INCLUDE_MODELVARIABLE_H_

#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>

template <int dim>
class modelVariable
{
public:
  dealii::VectorizedArray<double>                         scalarValue;
  dealii::Tensor<1, dim, dealii::VectorizedArray<double>> scalarGrad;
  dealii::Tensor<2, dim, dealii::VectorizedArray<double>> scalarHess;

  dealii::Tensor<1, dim, dealii::VectorizedArray<double>> vectorValue;
  dealii::Tensor<2, dim, dealii::VectorizedArray<double>> vectorGrad;
  dealii::Tensor<3, dim, dealii::VectorizedArray<double>> vectorHess;
};

template <int dim>
class modelResidual
{
public:
  dealii::VectorizedArray<double>                         scalarValueResidual;
  dealii::Tensor<1, dim, dealii::VectorizedArray<double>> scalarGradResidual;

  dealii::Tensor<1, dim, dealii::VectorizedArray<double>> vectorValueResidual;
  dealii::Tensor<2, dim, dealii::VectorizedArray<double>> vectorGradResidual;
};

struct variable_info
{
  bool         is_scalar;
  unsigned int scalar_or_vector_index;
  unsigned int global_var_index;
  bool         need_value;
  bool         need_gradient;
  bool         need_hessian;
  bool         value_residual;
  bool         gradient_residual;
  bool         var_needed;
};

#endif /* INCLUDE_MODELVARIABLE_H_ */
