
#ifndef PBasisSet_HH
#define PBasisSet_HH

#include "./PBasisSetBase.hh"
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
  class PBasisSet
  {
  public:
    std::string
    name() const
    {
      return (*ptr).name();
    }

    std::string
    description() const
    {
      return (*ptr)._description();
    }

    int
    size() const
    {
      return (*ptr).size();
    }

    virtual void
    resize(int N)
    {
      (*ptr).resize(N);
    }

    virtual int
    max_size()
    {
      // default to (essentially) no limit
      return (*ptr).max_size();
    }

    // Output PFunction for individual basis functions:

    PFunction<InType, OutType>
    basis_function(int term)
    {
      return (*ptr).basis_function(term);
    };

    // ----------------------------------------------------------
    // Use these functions if you want to evaluate a single value
    OutType
    operator()(int term, const InType &var)
    {
      return (*ptr)(term, var);
    }

    OutType
    grad(int term, const InType &var)
    {
      return (*ptr).grad(term, var);
    }

    OutType
    hess(int term, const InType &var)
    {
      return (*ptr).hess(term, var);
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
      return (*ptr).eval(var);
    }

    virtual const std::vector<OutType> &
    eval_grad(const InType &var)
    {
      return (*ptr).eval_grad(var);
    }

    virtual const std::vector<OutType> &
    eval_hess(const InType &var)
    {
      return (*ptr).eval_hess(var);
    }

    // Getters for individual terms
    OutType
    operator()(int term) const
    {
      return (*ptr)(term);
    }

    OutType
    grad(int term) const
    {
      return (*ptr).grad(term);
    }

    OutType
    hess(int term) const
    {
      return (*ptr).hess(term);
    }

    // Getters returning vector containing all terms
    const std::vector<OutType> &
    operator()() const
    {
      return (*ptr)();
    }

    const std::vector<OutType> &
    grad() const
    {
      return (*ptr).grad();
    }

    const std::vector<OutType> &
    hess() const
    {
      return (*ptr).hess();
    }

    // PFunction unique members ------------------------------------------

    PBasisSet &
    operator=(const PBasisSet &RHS)
    {
      if (ptr != NULL)
        delete ptr;
      ptr = RHS.ptr->clone();
      return *this;
    }

    template <class T>
    PBasisSet &
    operator=(const T &RHS)
    {
      RHS.is_derived_from_PBasisSetBase();

      if (ptr != NULL)
        delete ptr;
      ptr = RHS.clone();
      return *this;
    }

    // If you use this, PBasisSet becomes the 'owner' of the function RHS points to
    //    and it will delete it
    PBasisSet &
    set(PBasisSet *RHS)
    {
      if (RHS == NULL)
        {
          std::cout << "Error in PBasisSet::set. RHS == NULL." << std::endl;
          exit(1);
        }
      if (ptr != NULL)
        delete ptr;
      ptr = RHS;
      return *this;
    }

    PBasisSet()
    {
      ptr = NULL;
    }

    PBasisSet(const PBasisSet &RHS)
    {
      if (RHS.ptr != NULL)
        ptr = RHS.ptr->clone();
      else
        ptr = NULL;
    }

    template <class T>
    PBasisSet(const T &RHS)
    {
      RHS.is_derived_from_PBasisSetBase();

      ptr = RHS.clone();
    }

    ~PBasisSet()
    {
      if (ptr != NULL)
        delete ptr;
    }

  private:
    PBasisSetBase<InType, OutType> *ptr;
  };

} // namespace PRISMS

#endif