#ifndef INCLUDE_EQUATIONDEPENDECYPARSER_H_
#define INCLUDE_EQUATIONDEPENDECYPARSER_H_

# #include "dealIIheaders.h"


/**
* This is a class to parse the strings the user puts into the variableAttributeLoader to specify which variable values,
* gradients, hessians, etc are needed for each governing equation.
*/
class EquationDependencyParser
{
public:

  /**
  * Constructor, which currently is the only iterface to the outside world.
  */
  EquationDependencyParser(
      std::vector<std::string> sorted_dependencies_value_RHS,
      std::vector<std::string> sorted_dependencies_gradient_RHS,
      std::vector<std::string> sorted_dependencies_value_LHS,
      std::vector<std::string> sorted_dependencies_gradient_LHS);


    // All of the vectors of flags for what is needed
    std::vector<bool> need_value_explicit_RHS,
    need_gradient_explicit_RHS,
    need_hessian_explicit_RHS,
    need_value_nonexplicit_RHS,
    need_gradient_nonexplicit_RHS,
    need_hessian_nonexplicit_RHS,
    need_value_nonexplicit_LHS,
    need_gradient_nonexplicit_LHS,
    need_hessian_nonexplicit_LHS,
    need_value_change_nonexplicit_LHS,
    need_gradient_change_nonexplicit_LHS,
    need_hessian_change_nonexplicit_LHS;
};

#endif
