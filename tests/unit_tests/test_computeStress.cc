// Unit test(s) for the method "computeStress"

#include <iostream>
#include "test_computeStress.h"
#include "../../src/models/mechanics/computeStress.h"
//using namespace dealii;

#define problemDIM 3

//template <int dim>
void unitTest::test_computeStress(){
	bool pass = false;
	std::cout << "Testing 'computeStress...'" << std::endl;

	// Initialize CIJ
	dealii::Table<2, double> CIJ(2*problemDIM-1+problemDIM/3,2*problemDIM-1+problemDIM/3);
	CIJ.fill(1.0);

	if (problemDIM == 1) {
		CIJ[0][0] = 2.5;
	}
	else if (problemDIM == 2) {
		CIJ[0][0] = 6.8;
		CIJ[1][0] = 2.5;
		CIJ[2][0] = 4.0;
		CIJ[0][1] = 1.2;
		CIJ[1][1] = 10.1;
		CIJ[2][1] = 3.7;
		CIJ[0][2] = 9.3;
		CIJ[1][2] = 2.7;
		CIJ[2][2] = 8.8;
	}
	else if (problemDIM == 3) {
		CIJ[0][0] = 1.1;
		CIJ[1][1] = 7.7;
		CIJ[2][2] = 6.6;
		CIJ[3][3] = 3.3;
		CIJ[4][4] = 11.6;
		CIJ[5][5] = 19.5;
		CIJ[0][1]=CIJ[1][0] = 9.5;
		CIJ[0][2]=CIJ[2][0] = 2.1;
		CIJ[0][3]=CIJ[3][0] = 1.5;
		CIJ[0][4]=CIJ[4][0] = 9.2;
		CIJ[0][5]=CIJ[5][0] = 18.6;
		CIJ[1][2]=CIJ[2][1] = 5.6;
		CIJ[1][3]=CIJ[3][1] = 4.7;
		CIJ[1][4]=CIJ[4][1] = 6.4;
		CIJ[1][5]=CIJ[5][1] = 5.9;
		CIJ[2][3]=CIJ[3][2] = 15.5;
		CIJ[2][4]=CIJ[4][2] = 63.1;
		CIJ[2][5]=CIJ[5][2] = 50.0;
		CIJ[3][4]=CIJ[4][3] = 92.5;
		CIJ[3][5]=CIJ[5][3] = 1.3;
		CIJ[4][5]=CIJ[5][4] = 23.2;
	}

	// Initialize ux, a dim by dim 2nd rank tensor
	vectorgradType ux;
	if (problemDIM == 1) {
		ux[0][0] = 1.0;
	}
	else if (problemDIM == 2) {
		ux[0][0] = 1.0;
		ux[1][0] = 2.0;
		ux[0][1] = 3.0;
		ux[1][1] = 4.0;
	}
	else if (problemDIM == 3) {
		ux[0][0] = 1.0;
		ux[1][0] = 2.0;
		ux[2][0] = 3.0;
		ux[0][1] = 4.0;
		ux[1][1] = 5.0;
		ux[2][1] = 6.0;
		ux[0][2] = 7.0;
		ux[1][2] = 8.0;
		ux[2][2] = 9.0;
	}
	else {
			std::cerr << "Error: Number of dimensions must be between 1 and 3" << std::endl;
	}

	vectorgradType R;
	computeStress<problemDIM>(CIJ,ux,R);

	if (problemDIM == 1){
		if (R[0][0][0] - 2.5 < 1.0e-10) {pass = true;}
	}
	else if (problemDIM == 2){
		int pass_counter = 0;
		if (abs(R[0][0][0] - 58.1) < 1.0e-10) {pass_counter++;}
		if (abs(R[1][1][0] - 56.4) < 1.0e-10) {pass_counter++;}
		if (abs(R[0][1][0] - 62.8) < 1.0e-10) {pass_counter++;}
		if (abs(R[1][0][0] - 62.8) < 1.0e-10) {pass_counter++;}
		if (pass_counter == 4){pass = true;}
	}
	else if (problemDIM == 3){
		int pass_counter = 0;
		if (abs(R[0][0][0] - 292.1) < 1.0e-10) {pass_counter++;}
		if (abs(R[1][1][0] - 263.6) < 1.0e-10) {pass_counter++;}
		if (abs(R[2][2][0] - 1237.5) < 1.0e-10) {pass_counter++;}
		if (abs(R[1][2][0] - 1143.5) < 1.0e-10) {pass_counter++;}
		if (abs(R[2][1][0] - 1143.5) < 1.0e-10) {pass_counter++;}
		if (abs(R[0][2][0] - 2159.3) < 1.0e-10) {pass_counter++;}
		if (abs(R[2][0][0] - 2159.3) < 1.0e-10) {pass_counter++;}
		if (abs(R[0][1][0] - 865.3) < 1.0e-10) {pass_counter++;}
		if (abs(R[1][0][0] - 865.3) < 1.0e-10) {pass_counter++;}
		if (pass_counter == 9){pass = true;}
	}

	std::cout << "Test result for computeStress: " << pass << std::endl;

	std::cout << "Testing completed for 'computeStress'" << std::endl;
}




