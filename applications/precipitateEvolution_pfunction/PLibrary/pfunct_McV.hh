// created: 2017-1-11 12:28:14
// version: master
// url: https://github.com/bpuchala/IntegrationToolsWriter.git
// commit: 13e063c3ac8e8911a726a243fdbd68f291cc58cc

#ifndef pfunct_McV_HH
#define pfunct_McV_HH

#include "../../../include/IntegrationTools/PFunction.hh"
#include <cmath>
#include <cstdlib>

namespace PRISMS
{
  template <class VarContainer>
  class pfunct_McV_f : public PSimpleBase<VarContainer, double>
  {
    double
    eval([[maybe_unused]] const VarContainer &var) const override
    {
      return 1.0000000000000000e+00;
    }

  public:
    pfunct_McV_f()
    {
      this->_name = "pfunct_McV_f";
    }

    std::string
    csrc() const override
    {
      return "1.0000000000000000e+00";
    }

    std::string
    sym() const override
    {
      return "1.0";
    }

    std::string
    latex() const override
    {
      return "1.0";
    }

    pfunct_McV_f *
    clone() const override
    {
      return new pfunct_McV_f(*this);
    }
  };

  template <class VarContainer>
  class pfunct_McV : public PFuncBase<VarContainer, double>
  {
  public:
    typedef typename PFuncBase<VarContainer, double>::size_type size_type;

    PSimpleBase<VarContainer, double>   *_val;
    PSimpleBase<VarContainer, double>  **_grad_val;
    PSimpleBase<VarContainer, double> ***_hess_val;

    pfunct_McV()
    {
      construct();
    }

    pfunct_McV(const pfunct_McV &RHS)
      : PFuncBase<double *, double>(RHS)
    {
      construct(false);

      _val = RHS._val->clone();
    }

    pfunct_McV &
    operator=(pfunct_McV RHS)
    {
      using std::swap;

      swap(_val, RHS._val);

      return *this;
    }

    ~pfunct_McV()
    {
      delete _val;
    }

    pfunct_McV<VarContainer> *
    clone() const override
    {
      return new pfunct_McV<VarContainer>(*this);
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
      this->_name = "pfunct_McV";
      this->_var_name.clear();
      this->_var_name.push_back("c");
      this->_var_description.clear();
      this->_var_description.push_back("concentration");

      if (!allocate)
        return;

      _val = new pfunct_McV_f<VarContainer>();
    }
  };

} // namespace PRISMS
#endif
