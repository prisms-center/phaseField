template <int dim, typename T>
bool
unitTest<dim, T>::test_EquationDependencyParser_variables_and_residuals_needed()
{
  char buffer[100];

  std::cout << "\nTesting 'EquationDependencyParser'... " << std::endl;
  bool pass = true;

  // Declare all of the variables that go into the function call
  std::vector<std::string> sorted_dependencies_value_RHS,
    sorted_dependencies_gradient_RHS, sorted_dependencies_value_LHS,
    sorted_dependencies_gradient_LHS;

  // Populate the variable names and equation types
  std::vector<std::string> var_name;
  var_name.push_back("a");
  var_name.push_back("bc");

  std::vector<PDEType> var_eq_type;
  var_eq_type.push_back(EXPLICIT_TIME_DEPENDENT);
  var_eq_type.push_back(TIME_INDEPENDENT);

  // Populate the dependency strings
  sorted_dependencies_value_RHS.push_back("a, grad(a), bc");
  sorted_dependencies_value_RHS.push_back("hess(a),grad(bc)");

  sorted_dependencies_gradient_RHS.push_back("grad(bc), bc");
  sorted_dependencies_gradient_RHS.push_back("a,grad(bc)");

  sorted_dependencies_value_LHS.push_back("");
  sorted_dependencies_value_LHS.push_back("change(bc),grad(bc),grad(change(bc))");

  sorted_dependencies_gradient_LHS.push_back("");
  sorted_dependencies_gradient_LHS.push_back("bc");

  std::vector<bool> var_nonlinear;

  EquationDependencyParser equation_dependency_parser;

  equation_dependency_parser.parse(var_name,
                                   var_eq_type,
                                   sorted_dependencies_value_RHS,
                                   sorted_dependencies_gradient_RHS,
                                   sorted_dependencies_value_LHS,
                                   sorted_dependencies_gradient_LHS,
                                   var_nonlinear);

  // Check that the boolean vectors were generated correctly
  bool         result;
  unsigned int subtest_index = 0;

  // Check need_value_explicit_RHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_explicit_RHS.size() == 2)
    {
      if (equation_dependency_parser.eval_flags_explicit_RHS[0] &
            dealii::EvaluationFlags::values &&
          equation_dependency_parser.eval_flags_explicit_RHS[1] &
            dealii::EvaluationFlags::values)
        {
          result = true;
        }
    }
  std::cout << "Subtest " << subtest_index
            << " result for 'need_value_explicit_RHS': " << result << std::endl;
  pass = pass && result;

  // Check need_gradient_explicit_RHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_explicit_RHS.size() == 2)
    {
      if (equation_dependency_parser.eval_flags_explicit_RHS[0] &
            dealii::EvaluationFlags::gradients &&
          equation_dependency_parser.eval_flags_explicit_RHS[1] &
            dealii::EvaluationFlags::gradients)
        {
          result = true;
        }
    }

  std::cout << "Subtest " << subtest_index
            << " result for 'need_gradient_explicit_RHS': " << result << std::endl;
  pass = pass && result;

  // Check need_hessian_explicit_RHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_explicit_RHS.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_explicit_RHS[0] &
            dealii::EvaluationFlags::hessians) &&
          !(equation_dependency_parser.eval_flags_explicit_RHS[1] &
            dealii::EvaluationFlags::hessians))
        {
          result = true;
        }
    }

  std::cout << "Subtest " << subtest_index
            << " result for 'need_hessian_explicit_RHS': " << result << std::endl;
  pass = pass && result;

  // Check need_value_nonexplicit_RHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_nonexplicit_RHS.size() == 2)
    {
      if (equation_dependency_parser.eval_flags_nonexplicit_RHS[0] &
            dealii::EvaluationFlags::values &&
          !(equation_dependency_parser.eval_flags_nonexplicit_RHS[1] &
            dealii::EvaluationFlags::values))
        {
          result = true;
        }
    }

  std::cout << "Subtest " << subtest_index
            << " result for 'need_value_nonexplicit_RHS': " << result << std::endl;
  pass = pass && result;

  // Check need_gradient_nonexplicit_RHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_nonexplicit_RHS.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_nonexplicit_RHS[0] &
            dealii::EvaluationFlags::gradients) &&
          equation_dependency_parser.eval_flags_nonexplicit_RHS[1] &
            dealii::EvaluationFlags::gradients)
        {
          result = true;
        }
    }

  std::cout << "Subtest " << subtest_index
            << " result for 'need_gradient_nonexplicit_RHS': " << result << std::endl;
  pass = pass && result;

  // Check need_hessian_nonexplicit_RHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_nonexplicit_RHS.size() == 2)
    {
      if (equation_dependency_parser.eval_flags_nonexplicit_RHS[0] &
            dealii::EvaluationFlags::hessians &&
          !(equation_dependency_parser.eval_flags_nonexplicit_RHS[1] &
            dealii::EvaluationFlags::hessians))
        {
          result = true;
        }
    }

  std::cout << "Subtest " << subtest_index
            << " result for 'need_hessian_nonexplicit_RHS': " << result << std::endl;
  pass = pass && result;

  // Check need_value_nonexplicit_LHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_nonexplicit_LHS.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_nonexplicit_LHS[0] &
            dealii::EvaluationFlags::values) &&
          equation_dependency_parser.eval_flags_nonexplicit_LHS[1] &
            dealii::EvaluationFlags::values)
        {
          result = true;
        }
    }

  std::cout << "Subtest " << subtest_index
            << " result for 'need_value_nonexplicit_LHS': " << result << std::endl;
  pass = pass && result;

  // Check need_gradient_nonexplicit_LHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_nonexplicit_LHS.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_nonexplicit_LHS[0] &
            dealii::EvaluationFlags::gradients) &&
          equation_dependency_parser.eval_flags_nonexplicit_LHS[1] &
            dealii::EvaluationFlags::gradients)
        {
          result = true;
        }
    }

  std::cout << "Subtest " << subtest_index
            << " result for 'need_gradient_nonexplicit_LHS': " << result << std::endl;
  pass = pass && result;

  // Check need_hessian_nonexplicit_LHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_nonexplicit_LHS.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_nonexplicit_LHS[0] &
            dealii::EvaluationFlags::hessians) &&
          !(equation_dependency_parser.eval_flags_nonexplicit_LHS[1] &
            dealii::EvaluationFlags::hessians))
        {
          result = true;
        }
    }

  std::cout << "Subtest " << subtest_index
            << " result for 'need_hessian_nonexplicit_LHS': " << result << std::endl;
  pass = pass && result;

  // Check need_value_change_nonexplicit_LHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_change_nonexplicit_LHS.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_change_nonexplicit_LHS[0] &
            dealii::EvaluationFlags::values) &&
          equation_dependency_parser.eval_flags_change_nonexplicit_LHS[1] &
            dealii::EvaluationFlags::values)
        {
          result = true;
        }
    }

  std::cout << "Subtest " << subtest_index
            << " result for 'need_value_change_nonexplicit_LHS': " << result << std::endl;
  pass = pass && result;

  // Check need_gradient_change_nonexplicit_LHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_change_nonexplicit_LHS.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_change_nonexplicit_LHS[0] &
            dealii::EvaluationFlags::gradients) &&
          equation_dependency_parser.eval_flags_change_nonexplicit_LHS[1] &
            dealii::EvaluationFlags::gradients)
        {
          result = true;
        }
    }

  std::cout << "Subtest " << subtest_index
            << " result for 'need_gradient_change_nonexplicit_LHS': " << result
            << std::endl;
  pass = pass && result;

  // Check need_hessian_nonexplicit_LHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_change_nonexplicit_LHS.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_change_nonexplicit_LHS[0] &
            dealii::EvaluationFlags::hessians) &&
          !(equation_dependency_parser.eval_flags_change_nonexplicit_LHS[1] &
            dealii::EvaluationFlags::hessians))
        {
          result = true;
        }
    }
  std::cout << "Subtest " << subtest_index
            << " result for 'need_hessian_change_nonexplicit_LHS': " << result
            << std::endl;
  pass = pass && result;

  // Check need_value_residual_explicit_RHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_residual_explicit_RHS.size() == 2)
    {
      if (equation_dependency_parser.eval_flags_residual_explicit_RHS[0] &
            dealii::EvaluationFlags::values &&
          !(equation_dependency_parser.eval_flags_residual_explicit_RHS[1] &
            dealii::EvaluationFlags::values))
        {
          result = true;
        }
    }
  std::cout << "Subtest " << subtest_index
            << " result for 'need_value_residual_explicit_RHS': " << result << std::endl;
  pass = pass && result;

  // Check need_gradient_residual_explicit_RHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_residual_explicit_RHS.size() == 2)
    {
      if (equation_dependency_parser.eval_flags_residual_explicit_RHS[0] &
            dealii::EvaluationFlags::gradients &&
          !(equation_dependency_parser.eval_flags_residual_explicit_RHS[1] &
            dealii::EvaluationFlags::gradients))
        {
          result = true;
        }
    }
  std::cout << "Subtest " << subtest_index
            << " result for 'need_gradient_residual_explicit_RHS': " << result
            << std::endl;
  pass = pass && result;

  // Check need_value_residual_nonexplicit_RHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_residual_nonexplicit_RHS.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_residual_nonexplicit_RHS[0] &
            dealii::EvaluationFlags::values) &&
          equation_dependency_parser.eval_flags_residual_nonexplicit_RHS[1] &
            dealii::EvaluationFlags::values)
        {
          result = true;
        }
    }
  std::cout << "Subtest " << subtest_index
            << " result for 'need_value_residual_nonexplicit_RHS': " << result
            << std::endl;
  pass = pass && result;

  // Check need_gradient_residual_nonexplicit_RHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_residual_nonexplicit_RHS.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_residual_nonexplicit_RHS[0] &
            dealii::EvaluationFlags::gradients) &&
          equation_dependency_parser.eval_flags_residual_nonexplicit_RHS[1] &
            dealii::EvaluationFlags::gradients)
        {
          result = true;
        }
    }
  std::cout << "Subtest " << subtest_index
            << " result for 'need_gradient_residual_nonexplicit_RHS': " << result
            << std::endl;
  pass = pass && result;

  // Check need_value_residual_nonexplicit_LHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_residual_nonexplicit_LHS.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_residual_nonexplicit_LHS[0] &
            dealii::EvaluationFlags::values) &&
          equation_dependency_parser.eval_flags_residual_nonexplicit_LHS[1] &
            dealii::EvaluationFlags::values)
        {
          result = true;
        }
    }
  std::cout << "Subtest " << subtest_index
            << " result for 'need_value_residual_nonexplicit_LHS': " << result
            << std::endl;
  pass = pass && result;

  // Check need_gradient_residual_nonexplicit_LHS
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_residual_nonexplicit_LHS.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_residual_nonexplicit_LHS[0] &
            dealii::EvaluationFlags::gradients) &&
          equation_dependency_parser.eval_flags_residual_nonexplicit_LHS[1] &
            dealii::EvaluationFlags::gradients)
        {
          result = true;
        }
    }
  std::cout << "Subtest " << subtest_index
            << " result for 'need_gradient_residual_nonexplicit_LHS': " << result
            << std::endl;
  pass = pass && result;

  snprintf(buffer,
           sizeof(buffer),
           "Test result for 'EquationDependencyParser': %u\n",
           pass);
  std::cout << buffer;

  return pass;
}

template <int dim, typename T>
bool
unitTest<dim, T>::test_EquationDependencyParser_nonlinear()
{
  char buffer[100];

  std::cout << "\nTesting 'EquationDependencyParser'... " << std::endl;
  bool pass = true;

  // Declare all of the variables that go into the function call
  std::vector<std::string> sorted_dependencies_value_RHS,
    sorted_dependencies_gradient_RHS, sorted_dependencies_value_LHS,
    sorted_dependencies_gradient_LHS;

  // Populate the variable names and equation types
  std::vector<std::string> var_name;
  var_name.push_back("a");
  var_name.push_back("bc");
  var_name.push_back("def");
  var_name.push_back("ghij");
  var_name.push_back("klmno");

  std::vector<PDEType> var_eq_type;
  var_eq_type.push_back(EXPLICIT_TIME_DEPENDENT); // Not nonlinear by definition
  var_eq_type.push_back(TIME_INDEPENDENT);        // Crafted to be nonlinear due to RHS
  var_eq_type.push_back(
    TIME_INDEPENDENT); // Crafted to be nonlinear due to LHS (due to needed
                       // non-change version of governing variable)
  var_eq_type.push_back(TIME_INDEPENDENT); // Crafted to be nonlinear due to LHS (due to
                                           // needing a non-governing variable)
  var_eq_type.push_back(TIME_INDEPENDENT); // Crafted to be linear

  // Populate the dependency strings
  sorted_dependencies_value_RHS.push_back("bc, def, grad(bc)");
  sorted_dependencies_value_RHS.push_back("bc, def");
  sorted_dependencies_value_RHS.push_back("def");
  sorted_dependencies_value_RHS.push_back("a, ghij");
  sorted_dependencies_value_RHS.push_back("a, klmno");

  sorted_dependencies_gradient_RHS.push_back("hess(def), ghij");
  sorted_dependencies_gradient_RHS.push_back("grad(bc)");
  sorted_dependencies_gradient_RHS.push_back("");
  sorted_dependencies_gradient_RHS.push_back("grad(ghij)");
  sorted_dependencies_gradient_RHS.push_back("grad(klmno)");

  sorted_dependencies_value_LHS.push_back("");
  sorted_dependencies_value_LHS.push_back("change(bc)");
  sorted_dependencies_value_LHS.push_back("change(def)");
  sorted_dependencies_value_LHS.push_back("change(ghij), grad(def)");
  sorted_dependencies_value_LHS.push_back("change(klmno), grad(a)");

  sorted_dependencies_gradient_LHS.push_back("");
  sorted_dependencies_gradient_LHS.push_back("");
  sorted_dependencies_gradient_LHS.push_back("grad(def)");
  sorted_dependencies_gradient_LHS.push_back("");
  sorted_dependencies_gradient_LHS.push_back("hess(change(klmno))");

  std::vector<bool> var_nonlinear;

  EquationDependencyParser equation_dependency_parser;

  equation_dependency_parser.parse(var_name,
                                   var_eq_type,
                                   sorted_dependencies_value_RHS,
                                   sorted_dependencies_gradient_RHS,
                                   sorted_dependencies_value_LHS,
                                   sorted_dependencies_gradient_LHS,
                                   var_nonlinear);

  // Check that the boolean vectors were generated correctly
  bool         result;
  unsigned int subtest_index = 0;

  // Check explicit parabolic equation
  subtest_index++;
  result = false;
  if (var_nonlinear[0] == false)
    {
      result = true;
    }
  std::cout << "Subtest " << subtest_index
            << " result for an explicit parabolic equation: " << result << std::endl;
  pass = pass && result;

  // Check nonlinear time independent equation w/ nonlinearity due result for
  // nonlinearity due to needing a non-governing variable on the RHS
  subtest_index++;
  result = false;
  if (var_nonlinear[1] == true)
    {
      result = true;
    }
  std::cout << "Subtest " << subtest_index
            << " result for nonlinearity due to needing a non-governing "
               "variable on the RHS: "
            << result << std::endl;
  pass = pass && result;

  // Check nonlinear time independent equation w/ nonlinearity due to needed
  // non-change version of governing variable on the LHS
  subtest_index++;
  result = false;
  if (var_nonlinear[2] == true)
    {
      result = true;
    }
  std::cout << "Subtest " << subtest_index
            << " result for nonlinearity due to needed non-change version of "
               "governing variable on the LHS: "
            << result << std::endl;
  pass = pass && result;

  // Check nonlinear time independent equation w/ nonlinearity due to needing a
  // non-governing variable on the LHS
  subtest_index++;
  result = false;
  if (var_nonlinear[3] == true)
    {
      result = true;
    }
  std::cout << "Subtest " << subtest_index
            << " result for nonlinearity due to needing a non-governing "
               "variable on the LHS: "
            << result << std::endl;
  pass = pass && result;

  // Check linear time independent equation
  subtest_index++;
  result = false;
  if (var_nonlinear[4] == false)
    {
      result = true;
    }
  std::cout << "Subtest " << subtest_index
            << " result for linear time independent: " << result << std::endl;
  pass = pass && result;

  snprintf(buffer,
           sizeof(buffer),
           "Test result for 'EquationDependencyParser': %u\n",
           pass);
  std::cout << buffer;

  return pass;
}

template <int dim, typename T>
bool
unitTest<dim, T>::test_EquationDependencyParser_postprocessing()
{
  char buffer[100];

  std::cout << "\nTesting 'EquationDependencyParser'... " << std::endl;
  bool pass = true;

  // Declare all of the variables that go into the function call
  std::vector<std::string> sorted_dependencies_value, sorted_dependencies_gradient;

  // Populate the variable names and equation types
  std::vector<std::string> var_name;
  var_name.push_back("a");
  var_name.push_back("bc");

  std::vector<std::string> pp_var_name;
  pp_var_name.push_back("pp1");
  pp_var_name.push_back("pp2");
  pp_var_name.push_back("pp3");

  // Populate the dependency strings
  sorted_dependencies_value.push_back("bc, a, grad(bc)");
  sorted_dependencies_value.push_back("");
  sorted_dependencies_value.push_back("a");

  sorted_dependencies_gradient.push_back("hess(a), bc");
  sorted_dependencies_gradient.push_back("grad(bc)");
  sorted_dependencies_gradient.push_back("a");

  EquationDependencyParser equation_dependency_parser;

  equation_dependency_parser.pp_parse(var_name,
                                      pp_var_name,
                                      sorted_dependencies_value,
                                      sorted_dependencies_gradient);

  // Check that the boolean vectors were generated correctly
  bool         result;
  unsigned int subtest_index = 0;

  // Check pp_need_value
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_postprocess.size() == 2)
    {
      if (equation_dependency_parser.eval_flags_postprocess[0] &
            dealii::EvaluationFlags::values &&
          equation_dependency_parser.eval_flags_postprocess[1] &
            dealii::EvaluationFlags::values)
        {
          result = true;
        }
    }
  std::cout << "Subtest " << subtest_index << " result for 'pp_need_value': " << result
            << std::endl;
  pass = pass && result;

  // Check pp_need_gradient
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_postprocess.size() == 2)
    {
      if (!(equation_dependency_parser.eval_flags_postprocess[0] &
            dealii::EvaluationFlags::gradients) &&
          equation_dependency_parser.eval_flags_postprocess[1] &
            dealii::EvaluationFlags::gradients)
        {
          result = true;
        }
    }

  std::cout << "Subtest " << subtest_index << " result for 'pp_need_gradient': " << result
            << std::endl;
  pass = pass && result;

  // Check pp_need_hessian
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_postprocess.size() == 2)
    {
      if (equation_dependency_parser.eval_flags_postprocess[0] &
            dealii::EvaluationFlags::hessians &&
          !(equation_dependency_parser.eval_flags_postprocess[1] &
            dealii::EvaluationFlags::hessians))
        {
          result = true;
        }
    }

  std::cout << "Subtest " << subtest_index << " result for 'pp_need_hessian': " << result
            << std::endl;
  pass = pass && result;

  // Check pp_need_value_residual
  subtest_index++;
  result = false;
  if (equation_dependency_parser.eval_flags_residual_postprocess.size() == 3)
    {
      if (equation_dependency_parser.eval_flags_residual_postprocess[0] &
            dealii::EvaluationFlags::values &&
          !(equation_dependency_parser.eval_flags_residual_postprocess[1] &
            dealii::EvaluationFlags::values) &&
          equation_dependency_parser.eval_flags_residual_postprocess[2] &
            dealii::EvaluationFlags::values)
        {
          result = true;
        }
    }
  std::cout << "Subtest " << subtest_index
            << " result for 'pp_need_value_residual': " << result << std::endl;
  pass = pass && result;

  // Check pp_need_gradient_residual
  subtest_index++;
  result = false;

  if (equation_dependency_parser.eval_flags_residual_postprocess.size() == 3)
    {
      if (equation_dependency_parser.eval_flags_residual_postprocess[0] &
            dealii::EvaluationFlags::gradients &&
          equation_dependency_parser.eval_flags_residual_postprocess[1] &
            dealii::EvaluationFlags::gradients &&
          equation_dependency_parser.eval_flags_residual_postprocess[2] &
            dealii::EvaluationFlags::gradients)
        {
          result = true;
        }
    }
  std::cout << "Subtest " << subtest_index
            << " result for 'pp_need_gradient_residual': " << result << std::endl;
  pass = pass && result;

  snprintf(buffer,
           sizeof(buffer),
           "Test result for 'EquationDependencyParser': %u\n",
           pass);
  std::cout << buffer;

  return pass;
}
