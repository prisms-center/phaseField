
#ifndef PFuncBase_HH
#define PFuncBase_HH

#include "./PSimpleFunction.hh"
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
    std::string              _name;
    std::vector<std::string> _var_name;
    std::vector<std::string> _var_description;

  public:
    using size_type = std::vector<std::string>::size_type;

    PFuncBase()
    {}

    PFuncBase(const std::string              &name,
              const std::vector<std::string> &var_name,
              const std::vector<std::string> &var_description)
      : _name(name)
      , _var_name(var_name)
      , _var_description(var_description)
    {}

    virtual ~PFuncBase()
    {}

    std::string
    name()
    {
      return _name;
    }

    [[nodiscard]] size_type
    size() const
    {
      return _var_name.size();
    }

    std::vector<std::string>
    var_name()
    {
      return _var_name;
    }

    std::string
    var_name(size_type i)
    {
      return _var_name[i];
    }

    std::vector<std::string>
    var_description()
    {
      return _var_description;
    }

    std::string
    var_description(size_type i)
    {
      return _var_description[i];
    }

    void
    is_derived_from_PFuncBase() const
    {
      return;
    }

    virtual PFuncBase<VarContainer, OutType> *
    clone() const
    {
      return new PFuncBase<VarContainer, OutType>(*this);
    }

    virtual PSimpleFunction<VarContainer, OutType>
    simplefunction() const
    {
      undefined("PSimpleFunction<VarContainer, OutType> simplefunction() const");
      return PSimpleFunction<VarContainer, OutType>();
    }

    virtual PSimpleFunction<VarContainer, OutType>
    grad_simplefunction([[maybe_unused]] size_type di) const
    {
      undefined("PSimpleFunction<VarContainer, OutType> grad_simplefunction() const");
      return PSimpleFunction<VarContainer, OutType>();
    }

    virtual PSimpleFunction<VarContainer, OutType>
    hess_simplefunction([[maybe_unused]] size_type di,
                        [[maybe_unused]] size_type dj) const
    {
      undefined("PSimpleFunction<VarContainer, OutType> hess_simplefunction(size_type "
                "di, size_type dj) const");
      return PSimpleFunction<VarContainer, OutType>();
    }

    // ----------------------------------------------------------
    // Use these functions if you want to evaluate a single value
    virtual OutType
    operator()([[maybe_unused]] const VarContainer &var)
    {
      undefined("OutType operator()(const VarContainer &var)");
      return OutType();
    }

    virtual OutType
    grad([[maybe_unused]] const VarContainer &var, [[maybe_unused]] size_type di)
    {
      undefined("OutType grad(const VarContainer &var, size_type di)");
      return OutType();
    }

    virtual OutType
    hess([[maybe_unused]] const VarContainer &var,
         [[maybe_unused]] size_type           di,
         [[maybe_unused]] size_type           dj)
    {
      undefined("OutType hess(const VarContainer &var, size_type di, size_type dj)");
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

    virtual void
    eval_grad([[maybe_unused]] const VarContainer &var)
    {
      undefined("void eval_grad( const VarContainer &var)");
    }

    virtual void
    eval_hess([[maybe_unused]] const VarContainer &var)
    {
      undefined("void eval_hess( const VarContainer &var)");
    }

    virtual OutType
    operator()() const
    {
      undefined("OutType operator()");
      return OutType();
    }

    virtual OutType
    grad([[maybe_unused]] size_type di) const
    {
      undefined("OutType grad(size_type di)");
      return OutType();
    }

    virtual OutType
    hess([[maybe_unused]] size_type di, [[maybe_unused]] size_type dj) const
    {
      undefined("OutType hess(size_type di, size_type dj)");
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