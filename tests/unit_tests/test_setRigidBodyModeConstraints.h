// Unit test(s) for the method "setRigidBodyModeConstraints"
template <int dim, int degree>
class setRigidBodyModeConstraintsTest: public MatrixFreePDE<dim,degree>
{
	public:
	setRigidBodyModeConstraintsTest(userInputParameters _userInputs): MatrixFreePDE<dim,degree>(_userInputs) {
		//init the MatrixFreePDE class for testing
		this->initForTests();
	};

	void call_setRigidBodyModeConstraints(std::vector<int> rigidBodyModeComponents, unsigned int & num_constraints){

		this->setRigidBodyModeConstraints(rigidBodyModeComponents,this->constraintsOtherSet_nonconst[0],this->dofHandlersSet_nonconst[0]);

		// Calculate the number of constraints that were set
		num_constraints = this->constraintsOtherSet_nonconst[0]->n_constraints();
	};


	void setBCs(){};

	private:
	//RHS implementation for explicit solve
	void getRHS(const MatrixFree<dim,double> &data,
			std::vector<vectorType*> &dst,
			const std::vector<vectorType*> &src,
			const std::pair<unsigned int,unsigned int> &cell_range) const{};

	void residualRHS(const std::vector<modelVariable<dim> > & modelVarList,
			std::vector<modelResidual<dim> > & modelResidualsList,
			dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {};

	void residualLHS(const std::vector<modelVariable<dim> > & modelVarList,
			modelResidual<dim> & modelRes,
			dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {};

	void energyDensity(const std::vector<modelVariable<dim> > & modelVarList, const dealii::VectorizedArray<double> & JxW_value,
			dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) {};

};

template <int dim,typename T>
bool unitTest<dim,T>::test_setRigidBodyModeConstraints(std::vector<int> rigidBodyModeComponents, userInputParameters userInputs){

	bool pass = false;
	std::cout << "\nTesting 'setRigidBodyModeConstraints' with " << rigidBodyModeComponents.size() << " component(s) needing a constraint...'" << std::endl;

	//create test problem class object
	setRigidBodyModeConstraintsTest<dim,finiteElementDegree> test(userInputs);
	unsigned int num_constraints;
	test.call_setRigidBodyModeConstraints(rigidBodyModeComponents,num_constraints);

	// Add up the total number of constraints across all processors
	unsigned int global_num_constraints = rigidBodyModeComponents.size();
	Utilities::MPI::sum(global_num_constraints,MPI_COMM_WORLD);

	// Check if calculated value equals expected value
	if (global_num_constraints == rigidBodyModeComponents.size()){
		pass = true;
	}

	char buffer[100];
	sprintf (buffer, "Test result for 'setRigidBodyModeConstraints' with   %lu component(s) needing a constraint: %u\n", rigidBodyModeComponents.size(), pass);
	std::cout << buffer;

	return pass;
}
