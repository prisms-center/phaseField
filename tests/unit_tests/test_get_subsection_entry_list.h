// Unit test(s) for the method "parse_line"

template <int dim,typename T>
  bool unitTest<dim,T>::test_get_subsection_entry_list(){
  	bool pass = false;

    char buffer[100];

	std::cout << "\nTesting 'get_subsection_entry_list'... " << std::endl;

	//create test problem class object
    inputFileReader input_file_reader;
    std::vector<std::string> list_of_entries;
    list_of_entries = input_file_reader.get_subsection_entry_list("input_file_parser_test_file.txt","Equation","Variable type","SCALAR");

    if (list_of_entries.size() != 3){
        pass = false;
    }
    else {
        pass = (boost::iequals(list_of_entries[0],"SCALAR")
                    && boost::iequals(list_of_entries[1],"VECTOR")
                    && boost::iequals(list_of_entries[2],"SCALAR") );
    }

	sprintf (buffer, "Test result for 'get_subsection_entry_list': %u\n", pass);
	std::cout << buffer;

	return pass;
}
