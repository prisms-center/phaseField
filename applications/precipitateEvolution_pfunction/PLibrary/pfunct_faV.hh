// created: 2017-1-11 12:28:13
// version: master
// url: https://github.com/bpuchala/IntegrationToolsWriter.git
// commit: 13e063c3ac8e8911a726a243fdbd68f291cc58cc

#ifndef pfunct_faV_HH
#define pfunct_faV_HH

#include "../../../include/IntegrationTools/PFunction.hh"
#include <cmath>
#include <cstdlib>

namespace PRISMS
{
  template <class VarContainer>
  class pfunct_faV_f : public PSimpleBase<VarContainer, double>
  {
    double
    eval(const VarContainer &var) const override
    {
      return 5.1622000000000003e+00 * (var[0] * var[0]) +
             -2.7374999999999998e+00 * (var[0] * var[0] * var[0]) +
             -4.7759999999999998e+00 * var[0] +
             1.3687000000000000e+00 * ((var[0] * var[0]) * (var[0] * var[0])) -
             1.6704000000000001e+00;
    }

  public:
    pfunct_faV_f()
    {
      this->_name = "pfunct_faV_f";
    }

    std::string
    csrc() const override
    {
      return " 5.1622000000000003e+00*(var[0]*var[0])+-2.7374999999999998e+00*(var[0]*"
             "var[0]*var[0])+-4.7759999999999998e+00*var[0]+1.3687000000000000e+00*((var["
             "0]*var[0])*(var[0]*var[0]))-1.6704000000000001e+00";
    }

    std::string
    sym() const override
    {
      return "-1.6704+(1.3687)*c^4+(5.1622)*c^2-(2.7375)*c^3-(4.776)*c";
    }

    std::string
    latex() const override
    {
      return "-1.6704+{(5.1622)} c^{2}-{(4.776)} c-{(2.7375)} c^{3}+{(1.3687)} c^{4}";
    }

    pfunct_faV_f *
    clone() const override
    {
      return new pfunct_faV_f(*this);
    }
  };

  template <class VarContainer>
  class pfunct_faV_grad_0 : public PSimpleBase<VarContainer, double>
  {
    double
    eval(const VarContainer &var) const override
    {
      return 5.4748000000000001e+00 * (var[0] * var[0] * var[0]) +
             -8.2125000000000004e+00 * (var[0] * var[0]) +
             1.0324400000000001e+01 * var[0] - 4.7759999999999998e+00;
    }

  public:
    pfunct_faV_grad_0()
    {
      this->_name = "pfunct_faV_grad_0";
    }

    std::string
    csrc() const override
    {
      return " 5.4748000000000001e+00*(var[0]*var[0]*var[0])+-8.2125000000000004e+00*("
             "var[0]*var[0])+1.0324400000000001e+01*var[0]-4.7759999999999998e+00";
    }

    std::string
    sym() const override
    {
      return "-4.776+(10.3244)*c+(5.4748)*c^3-(8.2125)*c^2";
    }

    std::string
    latex() const override
    {
      return "-4.776+{(5.4748)} c^{3}-{(8.2125)} c^{2}+{(10.3244)} c";
    }

    pfunct_faV_grad_0 *
    clone() const override
    {
      return new pfunct_faV_grad_0(*this);
    }
  };

  template <class VarContainer>
  class pfunct_faV_hess_0_0 : public PSimpleBase<VarContainer, double>
  {
    double
    eval(const VarContainer &var) const override
    {
      return 1.6424399999999999e+01 * (var[0] * var[0]) +
             -1.6425000000000001e+01 * var[0] + 1.0324400000000001e+01;
    }

  public:
    pfunct_faV_hess_0_0()
    {
      this->_name = "pfunct_faV_hess_0_0";
    }

    std::string
    csrc() const override
    {
      return " 1.6424399999999999e+01*(var[0]*var[0])+-1.6425000000000001e+01*var[0]+1."
             "0324400000000001e+01";
    }

    std::string
    sym() const override
    {
      return "10.3244-(16.425)*c+(16.4244)*c^2";
    }

    std::string
    latex() const override
    {
      return "10.3244+{(16.4244)} c^{2}-{(16.425)} c";
    }

    pfunct_faV_hess_0_0 *
    clone() const override
    {
      return new pfunct_faV_hess_0_0(*this);
    }
  };

  template <class VarContainer>
  class pfunct_faV : public PFuncBase<VarContainer, double>
  {
  public:
    typedef typename PFuncBase<VarContainer, double>::size_type size_type;

    PSimpleBase<VarContainer, double>   *_val;
    PSimpleBase<VarContainer, double>  **_grad_val;
    PSimpleBase<VarContainer, double> ***_hess_val;

    pfunct_faV()
    {
      construct();
    }

    pfunct_faV(const pfunct_faV &RHS)
      : PFuncBase<double *, double>(RHS)
    {
      construct(false);

      _val            = RHS._val->clone();
      _grad_val[0]    = RHS._grad_val[0]->clone();
      _hess_val[0][0] = RHS._hess_val[0][0]->clone();
    }

    pfunct_faV &
    operator=(pfunct_faV RHS)
    {
      using std::swap;

      swap(_val, RHS._val);
      swap(_grad_val[0], RHS._grad_val[0]);
      swap(_hess_val[0][0], RHS._hess_val[0][0]);

      return *this;
    }

    ~pfunct_faV()
    {
      delete _val;

      delete _grad_val[0];
      delete[] _grad_val;

      delete _hess_val[0][0];
      delete[] _hess_val[0];
      delete[] _hess_val;
    }

    pfunct_faV<VarContainer> *
    clone() const override
    {
      return new pfunct_faV<VarContainer>(*this);
    }

    PSimpleFunction<VarContainer, double>
    simplefunction() const override
    {
      return PSimpleFunction<VarContainer, double>(*_val);
    }

    PSimpleFunction<VarContainer, double>
    grad_simplefunction(size_type di) const override
    {
      return PSimpleFunction<VarContainer, double>(*_grad_val[di]);
    }

    PSimpleFunction<VarContainer, double>
    hess_simplefunction(size_type di, size_type dj) const override
    {
      return PSimpleFunction<VarContainer, double>(*_hess_val[di][dj]);
    }

    double
    operator()(const VarContainer &var) override
    {
      return (*_val)(var);
    }

    double
    grad(const VarContainer &var, size_type di) override
    {
      return (*_grad_val[di])(var);
    }

    double
    hess(const VarContainer &var, size_type di, size_type dj) override
    {
      return (*_hess_val[di][dj])(var);
    }

    void
    eval(const VarContainer &var) override
    {
      (*_val)(var);
    }

    void
    eval_grad(const VarContainer &var) override
    {
      (*_grad_val[0])(var);
    }

    void
    eval_hess(const VarContainer &var) override
    {
      (*_hess_val[0][0])(var);
    }

    double
    operator()() const override
    {
      return (*_val)();
    }

    double
    grad(size_type di) const override
    {
      return (*_grad_val[di])();
    }

    double
    hess(size_type di, size_type dj) const override
    {
      return (*_hess_val[di][dj])();
    }

  private:
    void
    construct(bool allocate = true)
    {
      this->_name = "pfunct_faV";
      this->_var_name.clear();
      this->_var_name.push_back("c");
      this->_var_description.clear();
      this->_var_description.push_back("concentration");

      _grad_val = new PSimpleBase<VarContainer, double> *[1];

      _hess_val    = new PSimpleBase<VarContainer, double> **[1];
      _hess_val[0] = new PSimpleBase<VarContainer, double> *[1];

      if (!allocate)
        return;

      _val = new pfunct_faV_f<VarContainer>();

      _grad_val[0] = new pfunct_faV_grad_0<VarContainer>();

      _hess_val[0][0] = new pfunct_faV_hess_0_0<VarContainer>();
    }
  };

} // namespace PRISMS
#endif
