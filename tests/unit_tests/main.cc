// Calls the unit tests for the PRISMS-PF code
// Orignal author: Stephen DeWitt (stvdwtt@umich.edu)

#include <iostream>
#include "test_computeStress.h"

int main()
{
	std::cout << "Beginning unit tests..." << std::endl;

	bool pass_1D = false;
	bool pass_2D = false;
	bool pass_3D = false;


	unitTest<1> computeStress_tester_1D;
	pass_1D = computeStress_tester_1D.test_computeStress();

	unitTest<2> computeStress_tester_2D;
	pass_2D = computeStress_tester_2D.test_computeStress();

	unitTest<3> computeStress_tester_3D;
	pass_3D = computeStress_tester_3D.test_computeStress();

	int tests_passed = pass_1D+pass_2D+pass_3D;

	std::cout << std::endl;
	std::cout << "Number of tests passed: " << tests_passed << std::endl;

	return 0;
}





