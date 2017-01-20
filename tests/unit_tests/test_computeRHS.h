// Unit test(s) for the method "computeRHS"
#include "../../include/matrixFreePDE.h"

template <int dim,typename T>
bool unitTest<dim,T>::test_computeRHS(){

	bool pass = false;
	std::cout << "Testing 'computeRHS' in " << dim << " dimension(s)...'" << std::endl;

	MatrixFreePDE<dim> test_problem;
//	test_problem.fields.push_back(Field<dim>(SCALAR, PARABOLIC, "c"));
//	test_problem.fields.push_back(Field<dim>(SCALAR, PARABOLIC, "n1"));
//	test_problem.fields.push_back(Field<dim>(SCALAR, PARABOLIC, "n2"));
//	test_problem.fields.push_back(Field<dim>(SCALAR, PARABOLIC, "n3"));
//	test_problem.fields.push_back(Field<dim>(VECTOR,  ELLIPTIC, "u"));
	//test_problem.init ();

	td::cout << "The Test for 'computeRHS' is incomplete. Ignore result until the test is completed." << std::endl;

	std::cout << "Test result for 'computeRHS' in " << dim << " dimension(s): " << pass << std::endl;

	return pass;
}





