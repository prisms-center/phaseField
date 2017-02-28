// Model Variables Class

#ifndef INCLUDE_MODELVARIABLE_H_
#define INCLUDE_MODELVARIABLE_H_

template<int dim>
class modelVariable
{
 public:
	modelVariable();
	dealii::VectorizedArray<double> scalarValue;
	dealii::Tensor<1, dim, dealii::VectorizedArray<double> > scalarGrad;
	dealii::Tensor<2,dim,dealii::VectorizedArray<double> > scalarHess;

	dealii::Tensor<1, dim, dealii::VectorizedArray<double> > vectorValue;
	dealii::Tensor<2, dim, dealii::VectorizedArray<double> > vectorGrad;
	dealii::Tensor<3, dim, dealii::VectorizedArray<double> > vectorHess;
};

//constructor
template<int dim>
modelVariable<dim>::modelVariable()
{

}

template<int dim>
class modelResidual
{
 public:
	modelResidual();
	dealii::VectorizedArray<double> scalarValueResidual;
	dealii::Tensor<1, dim, dealii::VectorizedArray<double> > scalarGradResidual;

	dealii::Tensor<1, dim, dealii::VectorizedArray<double> > vectorValueResidual;
	dealii::Tensor<2, dim, dealii::VectorizedArray<double> > vectorGradResidual;

};

//constructor
template<int dim>
modelResidual<dim>::modelResidual()
{

}

struct variable_info
{
	bool is_scalar;
	unsigned int scalar_or_vector_index;
	unsigned int global_var_index;
};

#endif /* INCLUDE_MODELVARIABLE_H_ */





