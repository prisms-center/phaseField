// pFunction: An interface for PFunctions in PRISMS-PF
// This class is used to declare and load PFunctions from a pre-existing
// PLibrary that was created using PRISMS IntegrationTools. The class is a
// wrapper for the PFunctions in IntegrationTools so that they can be used for
// vectorized arrays. This class also calls the IntegrationTools checkout
// function in the constructor to reduce the number of steps the user needs to
// take. Currently this is only implemented for scalar functions. Vector
// functions can be treated component by component.

#include <deal.II/base/vectorization.h>

#include "../../include/IntegrationTools/extern/PLibraryExtern.hh"
#include <string>

namespace PFunctions
{
  using scalarvalueType = dealii::VectorizedArray<double>;

  class pFunction
  {
  public:
    // Constructor, wraps the IntegrationTools checkout function
    pFunction(const std::string &function_name);

    // Returns the value of the function for a given input variable
    scalarvalueType val(scalarvalueType);

    // Returns one of first derivatives of the function for a given input
    // variable
    scalarvalueType
    grad(scalarvalueType, unsigned int);

    // Returns one of the second derivatives of the function for a given input
    // variable
    scalarvalueType
    hess(scalarvalueType, unsigned int, unsigned int);

  private:
    PRISMS::PFunction<double *, double> fun;
  };

  pFunction::pFunction(const std::string &function_name)
  {
    PRISMS::PLibrary::checkout(function_name, fun);
  }

  scalarvalueType
  pFunction::val(scalarvalueType var)
  {
    scalarvalueType fun_val;
    for (unsigned i = 0; i < var.size(); i++)
      {
        fun_val[i] = fun(&var[i]);
      }
    return fun_val;
  }

  scalarvalueType
  pFunction::grad(scalarvalueType var, unsigned int dir)
  {
    scalarvalueType fun_grad;
    for (unsigned i = 0; i < var.size(); i++)
      {
        fun_grad[i] = fun.grad(&var[i], dir);
      }
    return fun_grad;
  }

  scalarvalueType
  pFunction::hess(scalarvalueType var, unsigned int dir1, unsigned int dir2)
  {
    scalarvalueType fun_hess;
    for (unsigned i = 0; i < var.size(); i++)
      {
        fun_hess[i] = fun.hess(&var[i], dir1, dir2);
      }
    return fun_hess;
  }

} // namespace PFunctions
