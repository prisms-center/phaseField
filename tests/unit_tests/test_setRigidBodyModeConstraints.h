// Unit test(s) for the method "setRigidBodyModeConstraints"
template <int dim>
class setRigidBodyModeConstraintsTest: public MatrixFreePDE<dim>
{
	public:
	setRigidBodyModeConstraintsTest(){
		//init the MatrixFreePDE class for testing
		this->initForTests();
	};

	void call_setRigidBodyModeConstraints(std::vector<int> rigidBodyModeComponents, unsigned int & num_constraints){

		this->setRigidBodyModeConstraints(rigidBodyModeComponents,this->constraintsOtherSet_nonconst[0],this->dofHandlersSet_nonconst[0]);

		// Calculate the number of constraints that were set
		num_constraints = this->constraintsOtherSet_nonconst[0]->n_constraints();
	};


 private:
	//RHS implementation for explicit solve
	  void getRHS(const MatrixFree<dim,double> &data,
		      std::vector<vectorType*> &dst,
		      const std::vector<vectorType*> &src,
		      const std::pair<unsigned int,unsigned int> &cell_range) const{};

};

template <int dim,typename T>
bool unitTest<dim,T>::test_setRigidBodyModeConstraints(std::vector<int> rigidBodyModeComponents){

	bool pass = false;
	std::cout << "\nTesting 'setRigidBodyModeConstraints' with " << rigidBodyModeComponents.size() << " component(s) needing a constraint...'" << std::endl;

	//create test problem class object
	setRigidBodyModeConstraintsTest<dim> test;
	unsigned int num_constraints;
	test.call_setRigidBodyModeConstraints(rigidBodyModeComponents,num_constraints);

	// Check if calculated value equals expected value
	if (num_constraints == rigidBodyModeComponents.size()){
		pass = true;
	}

	char buffer[100];
	sprintf (buffer, "Test result for 'setRigidBodyModeConstraints' with   %lu component(s) needing a constraint: %u\n", rigidBodyModeComponents.size(), pass);
	std::cout << buffer;

	return pass;
}
