// Unit test(s) for the method "getRHS"
#include "../../include/fields.h"

template <int dim,typename T>
bool unitTest<dim,T>::test_getRHS(){

	bool pass = false;
	std::cout << "Testing 'getRHS' in " << dim << " dimension(s)...'" << std::endl;

	// Initialize "valsScalar" and "valsVector"
	// valsScalar is of type std::map<std::string, typeScalar*>. typeScalar is of type dealii::FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,1,double>.
	// valsVector is of type std::map<std::string, typeVector*>. typeVector is of type dealii::FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,problemDIM,double>.
	std::map<std::string, typeScalar*>  valsScalar;
	std::map<std::string, typeVector*>  valsVector;
	std::map<std::string, unsigned int> valsIndex;

	// Start off with one scalar field, need to initialize a field
	//valsScalar[fields[0].name]=new typeScalar(data,0);
	std::vector<Field<dim> >                  fields;
	fields.push_back(Field<dim>(SCALAR, PARABOLIC, "c"));
	fields.push_back(Field<dim>(SCALAR, PARABOLIC, "n1"));
	fields.push_back(Field<dim>(SCALAR, PARABOLIC, "n2"));
	fields.push_back(Field<dim>(SCALAR, PARABOLIC, "n3"));
	fields.push_back(Field<dim>(VECTOR,  ELLIPTIC, "u"));

	// Initialize the MatrixFree object
	dealii::MatrixFree<dim,double> data;


	std::cout << "The Test for 'getRHS' is incomplete. Ignore result until the test is completed." << std::endl;

	std::cout << "Test result for 'getRHS' in " << dim << " dimension(s): " << pass << std::endl;

	return pass;
}
