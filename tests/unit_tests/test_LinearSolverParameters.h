// Unit tests for the class "nonlinearSolverParameters"
#include "../../include/SolverParameters.h"

template <int dim, typename T>
bool
unitTest<dim, T>::test_LinearSolverParameters()
{
  char buffer[100];

  std::cout << "\nTesting 'LinearSolverParameters'... " << std::endl;

  // Create a nonlinearSolverParameters object
  LinearSolverParameters test_object;

  test_object.loadParameters(2, ABSOLUTE_RESIDUAL, 1.0e-3, 123);
  test_object.loadParameters(5, RELATIVE_RESIDUAL_CHANGE, 1.0e-4, 124);

  // Subtests
  unsigned int subtest_index = 0;
  bool         result;
  bool         pass = true;

  // Subtest 1
  subtest_index++;
  result = false;
  if (test_object.getMaxIterations(2) == 123 && test_object.getMaxIterations(5) == 124)
    {
      result = true;
    }
  std::cout << "Subtest " << subtest_index << " result for 'getMaxIterations': " << result
            << std::endl;

  pass = pass && result;

  // Subtest 2
  subtest_index++;
  result = false;
  if (test_object.getToleranceType(2) == ABSOLUTE_RESIDUAL &&
      test_object.getToleranceType(5) == RELATIVE_RESIDUAL_CHANGE)
    {
      result = true;
    }
  std::cout << "Subtest " << subtest_index << " result for 'getToleranceType': " << result
            << std::endl;

  pass = pass && result;

  // Subtest 3
  subtest_index++;
  result = false;
  if (std::abs(test_object.getToleranceValue(2) - 1.0e-3) < 1.0e-12 &&
      std::abs(test_object.getToleranceValue(5) - 1.0e-4) < 1.0e-12)
    {
      result = true;
    }
  std::cout << "Subtest " << subtest_index
            << " result for 'getToleranceValue': " << result << std::endl;

  pass = pass && result;

  snprintf(buffer,
           sizeof(buffer),
           "Test result for 'LinearSolverParameters': %u\n",
           pass);
  std::cout << buffer;

  return pass;
}
