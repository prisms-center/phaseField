// Calls the unit tests for the PRISMS-PF code
// Orignal author: Stephen DeWitt (stvdwtt@umich.edu)

//#include "../../include/dealIIheaders.h"
//#include "../../src/models/coupled/coupledCHACMechanics.h"
#include <iostream>
#include "test_computeStress.h"
//using namespace dealii;


int main()
{
	std::cout << "Beginning unit tests..." << std::endl;

	const int dim = 1;
	unitTest computeStress_tester;

	//for (int dim = 1; dim < 4; dim++){

	computeStress_tester.test_computeStress();
	//}

	return 0;
}

