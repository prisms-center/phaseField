// Unit test(s) for the method "computeInvM"
template <int dim, int degree>
class testInvM : public MatrixFreePDE<dim, degree>
{
public:
  testInvM(userInputParameters<dim> _userInputs)
    : MatrixFreePDE<dim, degree>(_userInputs)
  {
    // Initialize the scalar test field object (needed for computeInvM())
    Field<dim> scalar_test_field(SCALAR, EXPLICIT_TIME_DEPENDENT, "c");
    this->fields.push_back(scalar_test_field);

    // Initialize the vector test field object (needed for computeInvM())
    Field<dim> vector_test_field(VECTOR, EXPLICIT_TIME_DEPENDENT, "c");
    this->fields.push_back(vector_test_field);

    // init the MatrixFreePDE class for testing
    this->initForTests(this->fields);

    // call computeInvM()
    this->computeInvM();
    invMNormscalar = this->invMscalar.l2_norm();
    invMNormvector = this->invMvector.l2_norm();
  };

  ~testInvM()
  {
    this->matrixFreeObject.clear();
  };

  double invMNormscalar;

  double invMNormvector;

  void
  setBCs() {};

  // Function to set the initial conditions (in ICs_and_BCs.h)
  void
  setInitialCondition(const dealii::Point<dim> &p,
                      const unsigned int        index,
                      double                   &scalar_IC,
                      dealii::Vector<double>   &vector_IC) override {};

  // Function to set the non-uniform Dirichlet boundary conditions (in
  // ICs_and_BCs.h)
  void
  setNonUniformDirichletBCs(const dealii::Point<dim> &p,
                            const unsigned int        index,
                            const unsigned int        direction,
                            const double              time,
                            double                   &scalar_BC,
                            dealii::Vector<double>   &vector_BC) override {};

private:
  // RHS implementation for explicit solve
  void
  getRHS(const MatrixFree<dim, double>               &data,
         std::vector<vectorType *>                   &dst,
         const std::vector<vectorType *>             &src,
         const std::pair<unsigned int, unsigned int> &cell_range) const {};

  // Function to set the RHS of the governing equations for explicit time
  // dependent equations (in equations.h)
  void
  explicitEquationRHS(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc,
    [[maybe_unused]] dealii::VectorizedArray<double> element_volume) const override {};

  // Function to set the RHS of the governing equations for all other equations
  // (in equations.h)
  void
  nonExplicitEquationRHS(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc,
    [[maybe_unused]] dealii::VectorizedArray<double> element_volume) const override {};

  // Function to set the LHS of the governing equations (in equations.h)
  void
  equationLHS(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    dealii::Point<dim, dealii::VectorizedArray<double>>              q_point_loc,
    [[maybe_unused]] dealii::VectorizedArray<double> element_volume) const override {};
};

template <int dim, typename T>
bool
unitTest<dim, T>::test_computeInvM(int                      argc,
                                   char                   **argv,
                                   userInputParameters<dim> userInputs)
{
  bool pass          = false;
  bool pass_subtest1 = false;
  bool pass_subtest2 = false;

  char buffer[100];
  std::cout << "\nTesting 'computeInvM' in " << dim << " dimension(s)...'" << std::endl;

  // create test problem class object
  // userInputParameters userInputs;
  // userInputs.loadUserInput();

  testInvM<dim, 1> test(userInputs);

  // Subtest 1 - check invMscalar norm
  if (std::abs(test.invMNormscalar - 1700.0) < 1.0e-10)
    {
      pass_subtest1 = true;
    }
  snprintf(buffer,
           sizeof(buffer),
           "Subtest 1 scalar result for 'computeInvM' in %u dimension(s): %u\n",
           dim,
           pass_subtest1);
  std::cout << buffer;

  // Subtest 2 - check invMscalar norm
  if (std::abs(test.invMNormvector - 2404.16305603426098) < 1.0e-10)
    {
      pass_subtest2 = true;
    }
  snprintf(buffer,
           sizeof(buffer),
           "Subtest 2 vector result for 'computeInvM' in %u dimension(s): %u\n",
           dim,
           pass_subtest2);
  std::cout << buffer;

  // Check if all subtests passed
  if (pass_subtest2 && pass_subtest1)
    {
      pass = true;
    }
  snprintf(buffer,
           sizeof(buffer),
           "Test result for 'computeInvM' in %u dimension(s): %u\n",
           dim,
           pass);
  std::cout << buffer;

  return pass;
}
