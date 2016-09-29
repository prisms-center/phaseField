// Interface for PFunctions in PRISMS-PF
// Used to declare and load PFunctions from a pre-existing PLibrary

namespace pFun{

class pFunction
{
public:
	pFunction(std::string function_name);

	scalarvalueType val(scalarvalueType);

	scalarvalueType grad(scalarvalueType, unsigned int);

	scalarvalueType hess(scalarvalueType, unsigned int, unsigned int);

private:
	PRISMS::PFunction<double*, double> fun;
};

pFunction::pFunction(std::string function_name){
	PRISMS::PLibrary::checkout(function_name, fun);
}

scalarvalueType pFunction::val(scalarvalueType var){
	scalarvalueType fun_val;
	for (unsigned i=0; i < var.n_array_elements; i++){
		fun_val[i] = fun(&var[i]);
	}
	return fun_val;

}


scalarvalueType pFunction::grad(scalarvalueType var,unsigned int dir){
	scalarvalueType fun_grad;
	for (unsigned i=0; i <var.n_array_elements; i++){
		fun_grad[i] = fun.grad(&var[i],dir);
	}
	return fun_grad;

}


scalarvalueType pFunction::hess(scalarvalueType var,unsigned int dir1, unsigned int dir2){
	scalarvalueType fun_hess;
	for (unsigned i=0; i < var.n_array_elements; i++){
		fun_hess[i] = fun.hess(&var[i],dir1,dir2);
	}
	return fun_hess;

}

}
