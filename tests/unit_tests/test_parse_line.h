// Unit test(s) for the method "parse_line"

template <int dim, typename T>
bool
unitTest<dim, T>::test_parse_line()
{
  bool pass          = false;
  bool pass_subtest1 = false;
  bool pass_subtest2 = false;
  bool pass_subtest3 = false;
  bool pass_subtest4 = false;
  bool pass_subtest5 = false;

  bool found_entry;
  char buffer[100];

  std::cout << "\nTesting 'parse_line'... " << std::endl;

  // create test problem class object
  // userInputParameters userInputs;
  // userInputs.loadUserInput();

  dealii::ParameterHandler parameter_handler;
  variableAttributeLoader  variable_attributes;
  inputFileReader          input_file_reader("parameters_test.prm", variable_attributes);

  // Subtest 1
  std::string line       = "set Test entry = 3";
  std::string keyword    = "set";
  std::string entry_name = "Test entry";
  std::string out_string;
  found_entry = input_file_reader.parse_line(line, keyword, entry_name, out_string, true);

  std::cout << out_string << " " << found_entry << std::endl;

  if (found_entry == true && boost::iequals(out_string, "3"))
    {
      pass_subtest1 = true;
    }
  snprintf(buffer,
           sizeof(buffer),
           "Subtest 1 result for 'parse_line': %u\n",
           pass_subtest1);
  std::cout << buffer;

  // Subtest 2
  line        = "set Test entry = 3";
  keyword     = "bork";
  entry_name  = "Test entry";
  found_entry = input_file_reader.parse_line(line, keyword, entry_name, out_string, true);

  if (found_entry == false)
    {
      pass_subtest2 = true;
    }
  snprintf(buffer,
           sizeof(buffer),
           "Subtest 2 result for 'parse_line': %u\n",
           pass_subtest2);
  std::cout << buffer;

  // Subtest 3
  line        = "set Test entry = 3";
  keyword     = "set";
  entry_name  = "Test entry wrong";
  found_entry = input_file_reader.parse_line(line, keyword, entry_name, out_string, true);

  if (found_entry == false)
    {
      pass_subtest3 = true;
    }
  snprintf(buffer,
           sizeof(buffer),
           "Subtest 3 result for 'parse_line': %u\n",
           pass_subtest3);
  std::cout << buffer;

  // Subtest 4
  line        = "set Test entry = 3";
  keyword     = "bork";
  entry_name  = "Test entry wrong";
  found_entry = input_file_reader.parse_line(line, keyword, entry_name, out_string, true);

  if (found_entry == false)
    {
      pass_subtest4 = true;
    }
  snprintf(buffer,
           sizeof(buffer),
           "Subtest 4 result for 'parse_line': %u\n",
           pass_subtest4);
  std::cout << buffer;

  // Subtest 5
  line        = "\tset Test entry = 3";
  keyword     = "set";
  entry_name  = "Test entry";
  found_entry = input_file_reader.parse_line(line, keyword, entry_name, out_string, true);

  if (found_entry == true && boost::iequals(out_string, "3"))
    {
      pass_subtest5 = true;
    }
  snprintf(buffer,
           sizeof(buffer),
           "Subtest 5 result for 'parse_line': %u\n",
           pass_subtest5);
  std::cout << buffer;

  // Check if all subtests passed
  if (pass_subtest4 && pass_subtest3 && pass_subtest2 && pass_subtest1)
    {
      pass = true;
    }

  snprintf(buffer, sizeof(buffer), "Test result for 'parse_line': %u\n", pass);
  std::cout << buffer;

  return pass;
}
