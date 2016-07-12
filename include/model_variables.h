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




