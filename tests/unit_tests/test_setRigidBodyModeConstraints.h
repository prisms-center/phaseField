// Unit test(s) for the method "setRigidBodyModeConstraints"
template <int dim, int degree>
class setRigidBodyModeConstraintsTest : public MatrixFreePDE<dim, degree>
{
public:
  setRigidBodyModeConstraintsTest(userInputParameters<dim> _userInputs)
    : MatrixFreePDE<dim, degree>(_userInputs)
  {
    // init the MatrixFreePDE class for testing
    Field<dim> scalar_test_field(SCALAR, EXPLICIT_TIME_DEPENDENT, "c");
    this->fields.push_back(scalar_test_field);
    this->initForTests(this->fields);
  };

  ~setRigidBodyModeConstraintsTest()
  {
    this->matrixFreeObject.clear();
  };

  void
  call_setRigidBodyModeConstraints(std::vector<int> rigidBodyModeComponents,
                                   unsigned int    &num_constraints)
  {
    this->setRigidBodyModeConstraints(rigidBodyModeComponents,
                                      this->constraintsOtherSet_nonconst[0],
                                      this->dofHandlersSet_nonconst[0]);

    // Calculate the number of constraints that were set
    num_constraints = this->constraintsOtherSet_nonconst[0]->n_constraints();
  };

  void
  setBCs() {};

  // Function to set the initial conditions (in ICs_and_BCs.h)
  void
  setInitialCondition(const dealii::Point<dim> &p,
                      const unsigned int        index,
                      double                   &scalar_IC,
                      dealii::Vector<double>   &vector_IC) {};

  // Function to set the non-uniform Dirichlet boundary conditions (in
  // ICs_and_BCs.h)
  void
  setNonUniformDirichletBCs(const dealii::Point<dim> &p,
                            const unsigned int        index,
                            const unsigned int        direction,
                            const double              time,
                            double                   &scalar_BC,
                            dealii::Vector<double>   &vector_BC) {};

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
    dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const {};

  // Function to set the RHS of the governing equations for all other equations
  // (in equations.h)
  void
  nonExplicitEquationRHS(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const {};

  // Function to set the LHS of the governing equations (in equations.h)
  void
  equationLHS(
    variableContainer<dim, degree, dealii::VectorizedArray<double>> &variable_list,
    dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc) const {};
};

template <int dim, typename T>
bool
unitTest<dim, T>::test_setRigidBodyModeConstraints(
  std::vector<int>         rigidBodyModeComponents,
  userInputParameters<dim> userInputs)
{
  bool pass = false;
  std::cout << "\nTesting 'setRigidBodyModeConstraints' with "
            << rigidBodyModeComponents.size() << " component(s) needing a constraint...'"
            << std::endl;

  // create test problem class object
  setRigidBodyModeConstraintsTest<dim, 1> test(userInputs);
  unsigned int                            num_constraints;
  test.call_setRigidBodyModeConstraints(rigidBodyModeComponents, num_constraints);

  // Add up the total number of constraints across all processors
  unsigned int global_num_constraints = rigidBodyModeComponents.size();
  Utilities::MPI::sum(global_num_constraints, MPI_COMM_WORLD);

  // Check if calculated value equals expected value
  if (global_num_constraints == rigidBodyModeComponents.size())
    {
      pass = true;
    }

  char buffer[100];
  snprintf(buffer,
           sizeof(buffer),
           "Test result for 'setRigidBodyModeConstraints' with   %lu "
           "component(s) needing a constraint: %u\n",
           rigidBodyModeComponents.size(),
           pass);
  std::cout << buffer;

  return pass;
}
