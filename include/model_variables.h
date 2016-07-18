// Model Variables Class

template<int dim>
class modelVariable
{
 public:
	modelVariable();
	scalarvalueType scalarValue;
	scalargradType scalarGrad;
	scalarhessType scalarHess;

	vectorvalueType vectorValue;
	vectorgradType vectorGrad;
	vectorhessType vectorHess;
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
	scalarvalueType scalarValueResidual;
	scalargradType scalarGradResidual;

	vectorvalueType vectorValueResidual;
	vectorgradType vectorGradResidual;

};

//constructor
template<int dim>
modelResidual<dim>::modelResidual()
{

}

template<int dim>
struct var_info
{
 //public:
	//var_info(MatrixFree<2,double> data, int index);
	//typeScalar varScalar;
	//typeVector varVector;
	bool is_scalar;
	unsigned int scalar_or_vector_index;
	unsigned int global_var_index;
	unsigned int global_field_index;


};







