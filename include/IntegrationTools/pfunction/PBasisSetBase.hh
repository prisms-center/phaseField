
#ifndef PBasisSetBase_HH
#define PBasisSetBase_HH

#include "./PFunction.hh"
#include <cstring>
#include <iostream>
#include <limits>
#include <vector>

namespace PRISMS
{
  /// Evaluate basis functions and their derivatives. Store and access the results.
  ///
  template <class InType, class OutType>
  class PBasisSetBase
  {
  public:
    typedef typename std::vector<OutType>::size_type size_type;

    std::string _name;
    std::string _description;

    std::vector<OutType> _val;
    std::vector<OutType> _grad_val;
    std::vector<OutType> _hess_val;

    std::string
    name() const
    {
      return _name;
    }

    std::string
    description() const
    {
      return _description;
    }

    size_type
    size() const
    {
      return _val.size();
    }

    PBasisSetBase(size_type N)
    {
      resize(N);
    }

    void
    is_derived_from_PBasisSetBase() const
    {
      return;
    }

    virtual void
    resize(size_type N)
    {
      _val.resize(N);
      _grad_val.resize(N);
      _hess_val.resize(N);
    }

    virtual ~PBasisSetBase()
    {}

    virtual size_type
    max_size() const
    {
      // default to (essentially) no limit
      return std::numeric_limits<size_type>::max();
    }

    virtual PBasisSetBase<InType, OutType> *
    clone() const
    {
      return new PBasisSetBase<InType, OutType>(*this);
    }

    virtual PFunction<InType, OutType>
    basis_function(size_type term) const
    {
      undefined("const PFunction<InType, OutType>& basis_function(size_type term) const");
      return PFunction<InType, OutType>();
    }

    // ----------------------------------------------------------
    // Use these functions if you want to evaluate a single value
    OutType
    operator()(size_type term, const InType &var)
    {
      return _val[term] = eval(term, var);
    }

    OutType
    grad(size_type term, const InType &var)
    {
      return _grad_val[term] = eval_grad(term, var);
    }

    OutType
    hess(size_type term, const InType &var)
    {
      return _hess_val[term] = eval_hess(term, var);
    }

    // ----------------------------------------------------------
    // Use these functions to evaluate several values, then use 'get' methods to access
    // results

    // By default, evaluate each individual term one-by-one.
    // These are virtual so derived classes may implement more efficient methods.
    //    Returns vector containing results
    virtual const std::vector<OutType> &
    eval(const InType &var)
    {
      for (size_type i = 0; i < _val.size(); i++)
        (*this)(i, var);
      return _val;
    }

    virtual const std::vector<OutType> &
    eval_grad(const InType &var)
    {
      for (size_type i = 0; i < _val.size(); i++)
        (*this).grad(i, var);
      return _grad_val;
    }

    virtual const std::vector<OutType> &
    eval_hess(const InType &var)
    {
      for (size_type i = 0; i < _val.size(); i++)
        (*this).hess(i, var);
      return _hess_val;
    }

    // Getters for individual terms
    OutType
    operator()(size_type term) const
    {
      return _val[term];
    }

    OutType
    grad(size_type term) const
    {
      return _grad_val[term];
    }

    OutType
    hess(size_type term) const
    {
      return _hess_val[term];
    }

    // Getters returning vector containing all terms
    const std::vector<OutType> &
    operator()() const
    {
      return _val;
    }

    const std::vector<OutType> &
    grad() const
    {
      return _grad_val;
    }

    const std::vector<OutType> &
    hess() const
    {
      return _hess_val;
    }

  private:
    /// ----------------------------------------------------------
    /// !!! Derived classes must define these functions !!!
    ///   Usually this is done using PBasisSetWriter
    virtual OutType
    eval(size_type term, const InType &var)
    {
      undefined("OutType PBasisSetBase::eval(size_type term, const InType &var)");
      return OutType();
    }

    virtual OutType
    eval_grad(size_type term, const InType &var)
    {
      undefined("OutType PBasisSetBase::eval_grad(size_type term, const InType &var)");
      return OutType();
    }

    virtual OutType
    eval_hess(size_type term, const InType &var)
    {
      undefined("OutType PBasisSetBase::eval_hess(size_type term, const InType &var)");
      return OutType();
    }

    void
    undefined(std::string fname) const
    {
      std::cout << "Error. In PBasisSetBase '" << _name << "'." << std::endl;
      std::cout << "   The member function '" << fname << "' has not been defined."
                << std::endl;
      exit(1);
    }
  };

} // namespace PRISMS

#endif