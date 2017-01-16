// created: 2017-1-11 12:28:13
// version: master
// url: https://github.com/bpuchala/IntegrationToolsWriter.git
// commit: 13e063c3ac8e8911a726a243fdbd68f291cc58cc

#ifndef pfunct_fbV_HH
#define pfunct_fbV_HH

#include <cmath>
#include <cstdlib>
#include "IntegrationTools/PFunction.hh"

namespace PRISMS
{
    template< class VarContainer>
    class pfunct_fbV_f : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  -5.9745999999999997e+00*var[0]+5.0000000000000000e+00*(var[0]*var[0])-1.5924000000000000e+00;
        }

    public:

        pfunct_fbV_f()
        {
            this->_name = "pfunct_fbV_f";
        }

        std::string csrc() const
        {
            return " -5.9745999999999997e+00*var[0]+5.0000000000000000e+00*(var[0]*var[0])-1.5924000000000000e+00";
        }

        std::string sym() const
        {
            return "-1.5924-(5.9746)*c+(5.0)*c^2";
        }

        std::string latex() const
        {
            return "-1.5924+{(5.0)} c^{2}-{(5.9746)} c";
        }

        pfunct_fbV_f* clone() const
        {
            return new pfunct_fbV_f(*this);
        }
    };

    template< class VarContainer>
    class pfunct_fbV_grad_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return  1.0000000000000000e+01*var[0]-5.9745999999999997e+00;
        }

    public:

        pfunct_fbV_grad_0()
        {
            this->_name = "pfunct_fbV_grad_0";
        }

        std::string csrc() const
        {
            return " 1.0000000000000000e+01*var[0]-5.9745999999999997e+00";
        }

        std::string sym() const
        {
            return "-5.9746+(10.0)*c";
        }

        std::string latex() const
        {
            return "-5.9746+{(10.0)} c";
        }

        pfunct_fbV_grad_0* clone() const
        {
            return new pfunct_fbV_grad_0(*this);
        }
    };

    template< class VarContainer>
    class pfunct_fbV_hess_0_0 : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 1.0000000000000000e+01;
        }

    public:

        pfunct_fbV_hess_0_0()
        {
            this->_name = "pfunct_fbV_hess_0_0";
        }

        std::string csrc() const
        {
            return "1.0000000000000000e+01";
        }

        std::string sym() const
        {
            return "10.0";
        }

        std::string latex() const
        {
            return "10.0";
        }

        pfunct_fbV_hess_0_0* clone() const
        {
            return new pfunct_fbV_hess_0_0(*this);
        }
    };

    template<class VarContainer>
    class pfunct_fbV : public PFuncBase< VarContainer, double>
    {
    public:
        
        typedef typename PFuncBase< VarContainer, double>::size_type size_type;

        PSimpleBase< VarContainer, double> *_val;
        PSimpleBase< VarContainer, double> **_grad_val;
        PSimpleBase< VarContainer, double> ***_hess_val;
        
        pfunct_fbV()
        {
            construct();
        }

        pfunct_fbV(const pfunct_fbV &RHS )
        {
            construct(false);
            
            _val = RHS._val->clone();
            _grad_val[0] = RHS._grad_val[0]->clone();
            _hess_val[0][0] = RHS._hess_val[0][0]->clone();
            
        }

        pfunct_fbV& operator=( pfunct_fbV RHS )
        {
            using std::swap;
            
            swap(_val, RHS._val);
            swap(_grad_val[0], RHS._grad_val[0]);
            swap(_hess_val[0][0], RHS._hess_val[0][0]);
            
            return *this;
        }

        ~pfunct_fbV()
        {
            delete _val;

            delete _grad_val[0];
            delete [] _grad_val;

            delete _hess_val[0][0];
            delete [] _hess_val[0];
            delete [] _hess_val;
        }

        pfunct_fbV<VarContainer>* clone() const
        {
            return new pfunct_fbV<VarContainer>(*this);
        }

        PSimpleFunction< VarContainer, double> simplefunction() const
        {
            return PSimpleFunction< VarContainer, double>( *_val );
        }

        PSimpleFunction< VarContainer, double> grad_simplefunction(size_type di) const
        {
            return PSimpleFunction< VarContainer, double>( *_grad_val[di] );
        }

        PSimpleFunction< VarContainer, double> hess_simplefunction(size_type di, size_type dj) const
        {
            return PSimpleFunction< VarContainer, double>( *_hess_val[di][dj] );
        }

        double operator()(const VarContainer &var)
        {
            return (*_val)(var);
        }

        double grad(const VarContainer &var, size_type di)
        {
            return (*_grad_val[di])(var);
        }

        double hess(const VarContainer &var, size_type di, size_type dj)
        {
            return (*_hess_val[di][dj])(var);
        }

        void eval(const VarContainer &var)
        {
            (*_val)(var);
        }

        void eval_grad(const VarContainer &var)
        {
            (*_grad_val[0])(var);
        }

        void eval_hess(const VarContainer &var)
        {
            (*_hess_val[0][0])(var);
        }

        double operator()() const
        {
            return (*_val)();
        }

        double grad(size_type di) const
        {
            return (*_grad_val[di])();
        }

        double hess(size_type di, size_type dj) const
        {
            return (*_hess_val[di][dj])();
        }

    private:
        void construct(bool allocate = true)
        {
            this->_name = "pfunct_fbV";
            this->_var_name.clear();
            this->_var_name.push_back("c");
            this->_var_description.clear();
            this->_var_description.push_back("concentration");
            
            _grad_val = new PSimpleBase< VarContainer, double>*[1];
            
            _hess_val = new PSimpleBase< VarContainer, double>**[1];
            _hess_val[0] = new PSimpleBase< VarContainer, double>*[1];
            
            if(!allocate) return;
            
            _val = new pfunct_fbV_f<VarContainer>();
            
            _grad_val[0] = new pfunct_fbV_grad_0<VarContainer>();
            
            _hess_val[0][0] = new pfunct_fbV_hess_0_0<VarContainer>();
        }

    };


}
#endif
