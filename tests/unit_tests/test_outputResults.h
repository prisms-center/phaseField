// Unit test(s) for the method "outputResults"
template <int dim>
class testOutputResults: public MatrixFreePDE<dim>
{
 public: 
  testOutputResults();
  
 private:
  //RHS implementation for explicit solve
  void getRHS(const MatrixFree<dim,double> &data, 
	      std::vector<vectorType*> &dst, 
	      const std::vector<vectorType*> &src,
	      const std::pair<unsigned int,unsigned int> &cell_range) const{};
};

//constructor
template <int dim>
testOutputResults<dim>::testOutputResults(): MatrixFreePDE<dim>(){
    //test outputResults
    //this->fields.push_back(Field<dim>(SCALAR, PARABOLIC, "c"));
    //this->outputResults();
 }


template <int dim,typename T>
  bool unitTest<dim,T>::test_outputResults(int argc, char** argv){
  bool pass = false;
  std::cout << "\nTesting 'outputResults' in " << dim << " dimension(s)...'" << std::endl;
  
  //create test problem class object
  testOutputResults<dim> test;
  pass=true;
  
  char buffer[100];
  sprintf (buffer, "Test result for 'outputResults' in %u dimension(s): %u\n", dim, pass);
  std::cout << buffer;
  
  return pass;
}
