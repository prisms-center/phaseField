template <int dim,typename T>
bool unitTest<dim,T>::test_EquationDependencyParser(){
    char buffer[100];

    std::cout << "\nTesting 'EquationDependencyParser'... " << std::endl;
    bool pass = true;


    // Declare all of the variables that go into the function call
    std::vector<std::string> sorted_dependencies_value_RHS, sorted_dependencies_gradient_RHS, sorted_dependencies_value_LHS, sorted_dependencies_gradient_LHS;

    // Populate the variable names and equation types
    std::vector<std::string> var_name;
    var_name.push_back("a");
    var_name.push_back("bc");

    std::vector<PDEType> var_eq_type;
    var_eq_type.push_back(PARABOLIC);
    var_eq_type.push_back(ELLIPTIC);


    // Populate the dependency strings
    sorted_dependencies_value_RHS.push_back("a, grad(a), bc");
    sorted_dependencies_value_RHS.push_back("hess(a),grad(bc)");

    sorted_dependencies_gradient_RHS.push_back("grad(bc), bc");
    sorted_dependencies_gradient_RHS.push_back("a,grad(bc)");

    sorted_dependencies_value_LHS.push_back("");
    sorted_dependencies_value_LHS.push_back("change(bc),grad(bc),grad(change(bc))");

    sorted_dependencies_gradient_LHS.push_back("");
    sorted_dependencies_gradient_LHS.push_back("change(a),bc");

    std::vector<bool> var_nonlinear;

    EquationDependencyParser equation_dependency_parser(
        var_name,
        var_eq_type,
        sorted_dependencies_value_RHS,
        sorted_dependencies_gradient_RHS,
        sorted_dependencies_value_LHS,
        sorted_dependencies_gradient_LHS,
        var_nonlinear);

    // Check that the boolean vectors were generated correctly
    bool result;
    unsigned int subtest_index = 0;

    // Check need_value_explicit_RHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_value_explicit_RHS.size() == 2){
        if (equation_dependency_parser.need_value_explicit_RHS[0] == true && equation_dependency_parser.need_value_explicit_RHS[1] == true){
            result = true;
        }
    }
    std::cout << "Subtest " << subtest_index << " result for 'need_value_explicit_RHS': " << result << std::endl;
    pass = pass && result;

    // Check need_gradient_explicit_RHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_value_explicit_RHS.size() == 2){
        if (equation_dependency_parser.need_gradient_explicit_RHS[0] == true && equation_dependency_parser.need_gradient_explicit_RHS[1] == true){
            result = true;
        }
    }

    std::cout << "Subtest " << subtest_index << " result for 'need_gradient_explicit_RHS': " << result << std::endl;
    pass = pass && result;

    // Check need_hessian_explicit_RHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_value_explicit_RHS.size() == 2){
        if (equation_dependency_parser.need_hessian_explicit_RHS[0] == false && equation_dependency_parser.need_hessian_explicit_RHS[1] == false){
            result = true;
        }
    }

    std::cout << "Subtest " << subtest_index << " result for 'need_hessian_explicit_RHS': " << result << std::endl;
    pass = pass && result;

    // Check need_value_nonexplicit_RHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_value_nonexplicit_RHS.size() == 2){
        if (equation_dependency_parser.need_value_nonexplicit_RHS[0] == true && equation_dependency_parser.need_value_nonexplicit_RHS[1] == false){
            result = true;
        }
    }

    std::cout << "Subtest " << subtest_index << " result for 'need_value_nonexplicit_RHS': " << result << std::endl;
    pass = pass && result;

    // Check need_gradient_nonexplicit_RHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_gradient_nonexplicit_RHS.size() == 2){
        if (equation_dependency_parser.need_gradient_nonexplicit_RHS[0] == false && equation_dependency_parser.need_gradient_nonexplicit_RHS[1] == true){
            result = true;
        }
    }

    std::cout << "Subtest " << subtest_index << " result for 'need_gradient_nonexplicit_RHS': " << result << std::endl;
    pass = pass && result;

    // Check need_hessian_nonexplicit_RHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_hessian_nonexplicit_RHS.size() == 2){
        if (equation_dependency_parser.need_hessian_nonexplicit_RHS[0] == true && equation_dependency_parser.need_hessian_nonexplicit_RHS[1] == false){
            result = true;
        }
    }

    std::cout << "Subtest " << subtest_index << " result for 'need_hessian_nonexplicit_RHS': " << result << std::endl;
    pass = pass && result;

    // Check need_value_nonexplicit_LHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_value_nonexplicit_LHS.size() == 2){
        if (equation_dependency_parser.need_value_nonexplicit_LHS[0] == false && equation_dependency_parser.need_value_nonexplicit_LHS[1] == true){
            result = true;
        }
    }

    std::cout << "Subtest " << subtest_index << " result for 'need_value_nonexplicit_LHS': " << result << std::endl;
    pass = pass && result;

    // Check need_gradient_nonexplicit_LHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_gradient_nonexplicit_LHS.size() == 2){
        if (equation_dependency_parser.need_gradient_nonexplicit_LHS[0] == false && equation_dependency_parser.need_gradient_nonexplicit_LHS[1] == true){
            result = true;
        }
    }

    std::cout << "Subtest " << subtest_index << " result for 'need_gradient_nonexplicit_LHS': " << result << std::endl;
    pass = pass && result;

    // Check need_hessian_nonexplicit_LHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_hessian_nonexplicit_LHS.size() == 2){
        if (equation_dependency_parser.need_hessian_nonexplicit_LHS[0] == false && equation_dependency_parser.need_hessian_nonexplicit_LHS[1] == false){
            result = true;
        }
    }

    std::cout << "Subtest " << subtest_index << " result for 'need_hessian_nonexplicit_LHS': " << result << std::endl;
    pass = pass && result;

    // Check need_value_change_nonexplicit_LHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_value_change_nonexplicit_LHS.size() == 2){
        if (equation_dependency_parser.need_value_change_nonexplicit_LHS[0] == true && equation_dependency_parser.need_value_change_nonexplicit_LHS[1] == true){
            result = true;
        }
    }

    std::cout << "Subtest " << subtest_index << " result for 'need_value_change_nonexplicit_LHS': " << result << std::endl;
    pass = pass && result;

    // Check need_gradient_change_nonexplicit_LHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_gradient_change_nonexplicit_LHS.size() == 2){
        if (equation_dependency_parser.need_gradient_change_nonexplicit_LHS[0] == false && equation_dependency_parser.need_gradient_change_nonexplicit_LHS[1] == true){
            result = true;
        }
    }

    std::cout << "Subtest " << subtest_index << " result for 'need_gradient_change_nonexplicit_LHS': " << result << std::endl;
    pass = pass && result;

    // Check need_hessian_nonexplicit_LHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_hessian_change_nonexplicit_LHS.size() == 2){
        if (equation_dependency_parser.need_hessian_change_nonexplicit_LHS[0] == false && equation_dependency_parser.need_hessian_change_nonexplicit_LHS[1] == false){
            result = true;
        }
    }
    std::cout << "Subtest " << subtest_index << " result for 'need_hessian_change_nonexplicit_LHS': " << result << std::endl;
    pass = pass && result;

    // Check need_value_residual_explicit_RHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_value_residual_explicit_RHS.size() == 2){
        if (equation_dependency_parser.need_value_residual_explicit_RHS[0] == true && equation_dependency_parser.need_value_residual_explicit_RHS[1] == false){
            result = true;
        }
    }
    std::cout << "Subtest " << subtest_index << " result for 'need_value_residual_explicit_RHS': " << result << std::endl;
    pass = pass && result;

    // Check need_gradient_residual_explicit_RHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_gradient_residual_explicit_RHS.size() == 2){
        if (equation_dependency_parser.need_gradient_residual_explicit_RHS[0] == true && equation_dependency_parser.need_gradient_residual_explicit_RHS[1] == false){
            result = true;
        }
    }
    std::cout << "Subtest " << subtest_index << " result for 'need_gradient_residual_explicit_RHS': " << result << std::endl;
    pass = pass && result;

    // Check need_value_residual_nonexplicit_RHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_value_residual_nonexplicit_RHS.size() == 2){
        if (equation_dependency_parser.need_value_residual_nonexplicit_RHS[0] == false && equation_dependency_parser.need_value_residual_nonexplicit_RHS[1] == true){
            result = true;
        }
    }
    std::cout << "Subtest " << subtest_index << " result for 'need_value_residual_nonexplicit_RHS': " << result << std::endl;
    pass = pass && result;

    // Check need_gradient_residual_nonexplicit_RHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_gradient_residual_nonexplicit_RHS.size() == 2){
        if (equation_dependency_parser.need_gradient_residual_nonexplicit_RHS[0] == false && equation_dependency_parser.need_gradient_residual_nonexplicit_RHS[1] == true){
            result = true;
        }
    }
    std::cout << "Subtest " << subtest_index << " result for 'need_gradient_residual_nonexplicit_RHS': " << result << std::endl;
    pass = pass && result;

    // Check need_value_residual_nonexplicit_LHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_value_residual_nonexplicit_LHS.size() == 2){
        if (equation_dependency_parser.need_value_residual_nonexplicit_LHS[0] == false && equation_dependency_parser.need_value_residual_nonexplicit_LHS[1] == true){
            result = true;
        }
    }
    std::cout << "Subtest " << subtest_index << " result for 'need_value_residual_nonexplicit_LHS': " << result << std::endl;
    pass = pass && result;

    // Check need_gradient_residual_nonexplicit_LHS
    subtest_index++;
    result = false;
    if (equation_dependency_parser.need_gradient_residual_nonexplicit_LHS.size() == 2){
        if (equation_dependency_parser.need_gradient_residual_nonexplicit_LHS[0] == false && equation_dependency_parser.need_gradient_residual_nonexplicit_LHS[1] == true){
            result = true;
        }
    }
    std::cout << "Subtest " << subtest_index << " result for 'need_gradient_residual_nonexplicit_LHS': " << result << std::endl;
    pass = pass && result;



    sprintf (buffer, "Test result for 'EquationDependencyParser': %u\n", pass);
	std::cout << buffer;

    return pass;

}
