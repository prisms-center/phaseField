// Unit test(s) for the method "computeInvM"
template <int dim>
class testInvM: public MatrixFreePDE<dim>
{
 public: 
  testInvM(){
    //init the MatrixFreePDE class for testing
    this->initForTests();

    //call computeInvM()
    this->computeInvM();
    invMNorm=this->invM.l2_norm();
  };
  double invMNorm;
  
 private:
  //RHS implementation for explicit solve
  void getRHS(const MatrixFree<dim,double> &data, 
	      std::vector<vectorType*> &dst, 
	      const std::vector<vectorType*> &src,
	      const std::pair<unsigned int,unsigned int> &cell_range) const{};

};

template <int dim,typename T>
  bool unitTest<dim,T>::test_computeInvM(int argc, char** argv){
  	bool pass = false;
	std::cout << "\nTesting 'computeInvM' in " << dim << " dimension(s)...'" << std::endl;
 
	//create test problem class object
	testInvM<dim> test;
	//check invM norm
	if ((test.invMNorm - 1700.0) < 1.0e-10) {pass=true;}
	char buffer[100];
	sprintf (buffer, "Test result for 'computeInvM' in   %u dimension(s): %u\n", dim, pass);
	std::cout << buffer;
	
	return pass;
}
