// Unit test(s) for the method "vectorLoad"
template <int dim, typename T>
bool unitTest<dim,T>::test_vectorLoad(T array[], int array_size, int num_array_elements){

	bool pass = true;
	std::cout << "Testing 'vectorLoad' for type '" << typeid(array).name() << "' ... " << std::endl;

	std::vector<T> test_vec;
	vectorLoad(array,array_size,test_vec);

	for (int i=0; i<num_array_elements; i++){
		if (test_vec[i] != array[i]){
			pass = false;
		}
	}

	std::cout << "Test result for 'vectorLoad' is " << pass << std::endl;

	return pass;
}

