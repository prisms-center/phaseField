// created: 2017-1-11 12:28:14
// version: master
// url: https://github.com/bpuchala/IntegrationToolsWriter.git
// commit: 13e063c3ac8e8911a726a243fdbd68f291cc58cc

#ifndef pfunct_Mn1V_HH
#define pfunct_Mn1V_HH

#include "../../../include/IntegrationTools/PFunction.hh"
#include <cmath>
#include <cstdlib>

namespace PRISMS
{
  template <class VarContainer>
  class pfunct_Mn1V_f : public PSimpleBase<VarContainer, double>
  {
    double
    eval([[maybe_unused]] const VarContainer &var) const override
    {
      return 1.0000000000000000e+02;
    }

  public:
    pfunct_Mn1V_f()
    {
      this->_name = "pfunct_Mn1V_f";
    }

    std::string
    csrc() const override
    {
      return "1.0000000000000000e+02";
    }

    std::string
    sym() const override
    {
      return "100.0";
    }

    std::string
    latex() const override
    {
      return "100.0";
    }

    pfunct_Mn1V_f *
    clone() const override
    {
      return new pfunct_Mn1V_f(*this);
    }
  };

  template <class VarContainer>
  class pfunct_Mn1V : public PFuncBase<VarContainer, double>
  {
  public:
    typedef typename PFuncBase<VarContainer, double>::size_type size_type;

    PSimpleBase<VarContainer, double>   *_val;
    PSimpleBase<VarContainer, double>  **_grad_val;
    PSimpleBase<VarContainer, double> ***_hess_val;

    pfunct_Mn1V()
    {
      construct();
    }

    pfunct_Mn1V(const pfunct_Mn1V &RHS)
      : PFuncBase<double *, double>(RHS)
    {
      construct(false);

      _val = RHS._val->clone();
    }

    pfunct_Mn1V &
    operator=(pfunct_Mn1V RHS)
    {
      using std::swap;

      swap(_val, RHS._val);

      return *this;
    }

    ~pfunct_Mn1V()
    {
      delete _val;
    }

    pfunct_Mn1V<VarContainer> *
    clone() const override
    {
      return new pfunct_Mn1V<VarContainer>(*this);
    }

    PSimpleFunction<VarContainer, double>
    simplefunction() const override
    {
      return PSimpleFunction<VarContainer, double>(*_val);
    }

    double
    operator()(const VarContainer &var) override
    {
      return (*_val)(var);
    }

    void
    eval(const VarContainer &var) override
    {
      (*_val)(var);
    }

    double
    operator()() const override
    {
      return (*_val)();
    }

  private:
    void
    construct(bool allocate = true)
    {
      this->_name = "pfunct_Mn1V";
      this->_var_name.clear();
      this->_var_name.push_back("n1");
      this->_var_description.clear();
      this->_var_description.push_back("concentration");

      if (!allocate)
        return;

      _val = new pfunct_Mn1V_f<VarContainer>();
    }
  };

} // namespace PRISMS
#endif
