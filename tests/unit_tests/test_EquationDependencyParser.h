template <int dim,typename T>
bool unitTest<dim,T>::test_EquationDependencyParser(){
    char buffer[100];

    std::cout << "\nTesting 'EquationDependencyParser'... " << std::endl;

    // Declare all of the variables that go into the function call
    std::vector<std::string> sorted_dependencies_value_RHS, sorted_dependencies_gradient_RHS, sorted_dependencies_value_LHS, sorted_dependencies_gradient_LHS;

    // Populate the variable names
    std::vector<std::string> var_name;
    var_name.push_back("a");
    var_name.push_back("bc");


    // Populate the dependency strings
    sorted_dependencies_value_RHS.push_back("a, grad(a), bc");
    sorted_dependencies_value_RHS.push_back("hess(a),grad(bc)");

    sorted_dependencies_gradient_RHS.push_back("a, grad(a), bc");
    sorted_dependencies_gradient_RHS.push_back("hess(a),grad(bc)");

    sorted_dependencies_value_LHS.push_back("a, change(a)");
    sorted_dependencies_value_LHS.push_back("hess(a),grad(bc)");

    sorted_dependencies_gradient_LHS.push_back("a, grad(a), bc");
    sorted_dependencies_gradient_LHS.push_back("hess(a),grad(bc)");


    parseEquationDependencies(
        sorted_dependencies_value_RHS,
        sorted_dependencies_gradient_RHS,
        sorted_dependencies_value_LHS,
        sorted_dependencies_gradient_LHS);


    std::cout << test.result << std::endl;

}
