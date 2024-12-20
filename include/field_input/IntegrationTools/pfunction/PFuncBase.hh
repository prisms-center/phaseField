#ifndef PFuncBase_HH
#define PFuncBase_HH

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace PRISMS
{
  /// A Base class for a function, including grad & hess
  ///
  template <class VarContainer, class OutType>
  class PFuncBase
  {
  protected:
    std::string _name;

  public:
    using size_type = std::vector<std::string>::size_type;

    PFuncBase() = default;

    PFuncBase(const std::string &name)
      : _name(name)
    {}

    virtual ~PFuncBase()
    {}

    std::string
    name()
    {
      return _name;
    }

    // ----------------------------------------------------------
    // Use these functions if you want to evaluate a single value
    virtual OutType
    operator()([[maybe_unused]] const VarContainer &var)
    {
      undefined("OutType operator()(const VarContainer &var)");
      return OutType();
    }

    // ----------------------------------------------------------
    // Use these functions to evaluate several values, then use 'get' methods to access
    // results
    virtual void
    eval([[maybe_unused]] const VarContainer &var)
    {
      undefined("void eval_grad( const VarContainer &var)");
    }

    virtual OutType
    operator()() const
    {
      undefined("OutType operator()");
      return OutType();
    }

  private:
    void
    undefined(std::string fname) const
    {
      std::string msg = "Error in PFuncBase '" + _name + "'.\n" +
                        "   The member function '" + fname + "' has not been defined.\n";
      throw std::runtime_error(msg);
    }
  };

} // namespace PRISMS

#endif