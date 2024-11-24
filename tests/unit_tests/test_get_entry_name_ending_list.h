// Unit test(s) for the method "get_entry_name_ending_list"

template <int dim, typename T>
bool
unitTest<dim, T>::test_get_entry_name_ending_list()
{
  bool pass = false;

  char buffer[100];

  std::cout << "\nTesting 'get_entry_name_ending_list'... " << std::endl;

  // create test problem class object
  variableAttributeLoader attribute_loader;
  const AttributesList    var_attributes = attribute_loader.get_var_attributes();
  const AttributesList    pp_attributes  = attribute_loader.get_pp_attributes();
  inputFileReader input_file_reader("parameters_test.prm", var_attributes, pp_attributes);
  std::vector<std::string> entry_name_ending_list;
  entry_name_ending_list =
    input_file_reader.get_entry_name_ending_list("get_entry_name_ending_list_tester.prm",
                                                 "set",
                                                 "Model constant");

  /*
  // Output the strings to see what, if anything, is going wrong
  for (unsigned int i=0; i < entry_name_ending_list.size(); i++){
      std::cout << entry_name_ending_list[i] << std::endl;
  }
  */

  if (entry_name_ending_list.size() != 5)
    {
      pass = false;
    }
  else
    {
      pass = (boost::iequals(entry_name_ending_list[0], "abc") &&
              boost::iequals(entry_name_ending_list[1], "dfg") &&
              boost::iequals(entry_name_ending_list[2], "i j k") &&
              boost::iequals(entry_name_ending_list[3], "qrs") &&
              boost::iequals(entry_name_ending_list[4], "t uv"));
    }

  snprintf(buffer,
           sizeof(buffer),
           "Test result for 'get_entry_name_ending_list': %u\n",
           pass);
  std::cout << buffer;

  return pass;
}
