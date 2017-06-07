// Calls the unit tests for the PRISMS-PF code
// Orignal author: Stephen DeWitt (stvdwtt@umich.edu)

#include <iostream>
#include "unitTest.h"
#include "../../include/userInputParameters.h"
#include "../../src/userInputParameters/loadUserInputs.cc"
#include "../../include/vectorBCFunction.h"
#include "../../src/utilities/vectorBCFunction.cc"
#include "../../include/initialConditions.h"
#include "initialConditions.cc"

int main(int argc, char **argv)
{
  //init MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, numbers::invalid_unsigned_int);

  // Load input
  userInputParameters<problemDIM> userInputs;
  userInputs.loadUserInput();

  std::cout << "Beginning unit tests..." << std::endl;
  bool pass = false;
  int tests_passed = 0, total_tests = 0;

  // Unit tests for the method "computeInvM"
  total_tests++;
  unitTest<2,double> computeInvM_tester_2D;
  pass = computeInvM_tester_2D.test_computeInvM(argc, argv, userInputs);
  tests_passed += pass;

  // Unit tests for the method "outputResults"
  total_tests++;
  unitTest<2,double> outputResults_tester_2D;
  pass = outputResults_tester_2D.test_outputResults(argc, argv, userInputs);
  tests_passed += pass;

  // Unit tests for the method "computeStress"
  total_tests++;
  unitTest<1,dealii::VectorizedArray<double>[1][1]> computeStress_tester_1D;
  pass = computeStress_tester_1D.test_computeStress();
  tests_passed += pass;

  total_tests++;
  unitTest<1,dealii::Table<2, double> > computeStress_tester_1DT;
  pass = computeStress_tester_1DT.test_computeStress();
  tests_passed += pass;

  total_tests++;
  unitTest<2,dealii::VectorizedArray<double>[3][3]> computeStress_tester_2D;
  pass = computeStress_tester_2D.test_computeStress();
  tests_passed += pass;

  total_tests++;
  unitTest<2,dealii::Table<2, double> > computeStress_tester_2DT;
  pass = computeStress_tester_2DT.test_computeStress();
  tests_passed += pass;

  total_tests++;
  unitTest<3,dealii::VectorizedArray<double>[6][6]> computeStress_tester_3D;
  pass = computeStress_tester_3D.test_computeStress();
  tests_passed += pass;

  total_tests++;
  unitTest<3,dealii::Table<2, double> > computeStress_tester_3DT;
  pass = computeStress_tester_3DT.test_computeStress();
  tests_passed += pass;

  // Unit tests for the method "getOutputTimeSteps" for all four types of spacing
  total_tests++;
  unitTest<2,double> getOutputTimeSteps_tester_eq;
  pass = getOutputTimeSteps_tester_eq.test_getOutputTimeSteps("EQUAL_SPACING",10,{},userInputs);
  tests_passed += pass;

  total_tests++;
  unitTest<2,double> getOutputTimeSteps_tester_log;
  pass = getOutputTimeSteps_tester_log.test_getOutputTimeSteps("LOG_SPACING",10,{}, userInputs);
  tests_passed += pass;

  total_tests++;
  unitTest<2,double> getOutputTimeSteps_tester_dec;
  pass = getOutputTimeSteps_tester_dec.test_getOutputTimeSteps("N_PER_DECADE",10,{}, userInputs);
  tests_passed += pass;

  total_tests++;
  unitTest<2,double> getOutputTimeSteps_tester_list;
  std::vector<unsigned int> userGivenTimeStepList = {0, 3, 55, 61};
  pass = getOutputTimeSteps_tester_list.test_getOutputTimeSteps("LIST",0,userGivenTimeStepList, userInputs);
  tests_passed += pass;

  // Unit tests for the method "setRigidBodyModeConstraints"
  total_tests++;
  unitTest<2,double> setRigidBodyModeConstraints_tester_null;
  std::vector<int> rigidBodyModeComponents;
  pass = setRigidBodyModeConstraints_tester_null.test_setRigidBodyModeConstraints(rigidBodyModeComponents, userInputs);
  tests_passed += pass;

  total_tests++;
  unitTest<2,double> setRigidBodyModeConstraints_tester_one;
  rigidBodyModeComponents.clear();
  rigidBodyModeComponents.push_back(0);
  pass = setRigidBodyModeConstraints_tester_one.test_setRigidBodyModeConstraints(rigidBodyModeComponents, userInputs);
  tests_passed += pass;

  total_tests++;
  unitTest<2,double> setRigidBodyModeConstraints_tester_three;
  rigidBodyModeComponents.clear();
  rigidBodyModeComponents.push_back(0);
  rigidBodyModeComponents.push_back(1);
  rigidBodyModeComponents.push_back(2);
  pass = setRigidBodyModeConstraints_tester_three.test_setRigidBodyModeConstraints(rigidBodyModeComponents, userInputs);
  tests_passed += pass;

  // Print out results
  char buffer[100];
  sprintf(buffer, "\n\nNumber of tests passed: %u/%u \n\n", tests_passed, total_tests);
  std::cout << buffer;

  // Write results to a file
  std::ofstream result_file ("unit_test_results.txt");
  result_file << tests_passed << std::endl;
  result_file << total_tests << std::endl;
  result_file.close();


  return 0;
}
