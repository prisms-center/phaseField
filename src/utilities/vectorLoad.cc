// =====================================================================
// LOADING INPUTS INTO INTERNAL VECTORS
// =====================================================================
// This could probably be simplified by using templates

void vectorLoad(bool in[], int array_size, std::vector<bool> & out){

	std::vector<bool> temp ( in, in+array_size / sizeof(bool) );
	out = temp;
}

void vectorLoad(double in[], int array_size, std::vector<double> & out){

	std::vector<double> temp ( in, in+array_size / sizeof(double) );
	out = temp;
}

void vectorLoad(int in[], int array_size, std::vector<int> & out){

	std::vector<int> temp ( in, in+array_size / sizeof(double) );
	out = temp;
}

void vectorLoad(unsigned int in[], int array_size, std::vector<unsigned int> & out){

	std::vector<unsigned int> temp ( in, in+array_size / sizeof(double) );
	out = temp;
}

void vectorLoad(std::string in[], int array_size, std::vector<std::string> & out){

	bool at_array_end = false;
	int size_loaded = 0;
	int index = 0;
	while (!at_array_end){
		if (size_loaded < array_size){
			out.push_back(in[index]);
			size_loaded += sizeof(in[index]);
			index++;
		}
		else {
			at_array_end = true;
		}
	}

}


