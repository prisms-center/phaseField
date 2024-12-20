#ifndef PSimpleBase_HH
#define PSimpleBase_HH

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

namespace PRISMS
{

  /// Base classes for functions that can be hard-coded,
  ///   then shared and used elsewhere

  /// A simple expression evaluator
  ///
  template <class VarContainer, class OutType>
  class PSimpleBase
  {
  protected:
    std::string _name;
    OutType     _val;

  public:
    virtual ~PSimpleBase() {};

    [[nodiscard]] std::string
    name() const
    {
      return _name;
    }

    OutType
    operator()(const VarContainer &var)
    {
      return _val = eval(var);
    }

    OutType
    operator()() const
    {
      return _val;
    }

  private:
    virtual OutType
    eval([[maybe_unused]] const VarContainer &var) const
    {
      undefined("OutType eval( const VarContainer &var)");
      return OutType();
    }

    void
    undefined(std::string fname) const
    {
      std::string msg = "Error in PSimpleBase '" + _name + "'.\n" +
                        "   The member function '" + fname + "' has not been defined.\n";
      throw std::runtime_error(msg);
    }
  };

} // namespace PRISMS

#endif