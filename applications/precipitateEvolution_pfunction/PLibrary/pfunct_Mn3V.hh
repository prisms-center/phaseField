// created: 2017-1-11 12:28:14
// version: master
// url: https://github.com/bpuchala/IntegrationToolsWriter.git
// commit: 13e063c3ac8e8911a726a243fdbd68f291cc58cc

#ifndef pfunct_Mn3V_HH
#define pfunct_Mn3V_HH

#include <cmath>
#include <cstdlib>
#include "IntegrationTools/PFunction.hh"

namespace PRISMS
{
    template< class VarContainer>
    class pfunct_Mn3V_f : public PSimpleBase< VarContainer, double>
    {
        double eval( const VarContainer &var) const
        {
            return 1.0000000000000000e+02;
        }

    public:

        pfunct_Mn3V_f()
        {
            this->_name = "pfunct_Mn3V_f";
        }

        std::string csrc() const
        {
            return "1.0000000000000000e+02";
        }

        std::string sym() const
        {
            return "100.0";
        }

        std::string latex() const
        {
            return "100.0";
        }

        pfunct_Mn3V_f* clone() const
        {
            return new pfunct_Mn3V_f(*this);
        }
    };

    template<class VarContainer>
    class pfunct_Mn3V : public PFuncBase< VarContainer, double>
    {
    public:
        
        typedef typename PFuncBase< VarContainer, double>::size_type size_type;

        PSimpleBase< VarContainer, double> *_val;
        PSimpleBase< VarContainer, double> **_grad_val;
        PSimpleBase< VarContainer, double> ***_hess_val;
        
        pfunct_Mn3V()
        {
            construct();
        }

        pfunct_Mn3V(const pfunct_Mn3V &RHS )
        {
            construct(false);
            
            _val = RHS._val->clone();
            
        }

        pfunct_Mn3V& operator=( pfunct_Mn3V RHS )
        {
            using std::swap;
            
            swap(_val, RHS._val);
            
            return *this;
        }

        ~pfunct_Mn3V()
        {
            delete _val;

        }

        pfunct_Mn3V<VarContainer>* clone() const
        {
            return new pfunct_Mn3V<VarContainer>(*this);
        }

        PSimpleFunction< VarContainer, double> simplefunction() const
        {
            return PSimpleFunction< VarContainer, double>( *_val );
        }

        double operator()(const VarContainer &var)
        {
            return (*_val)(var);
        }

        void eval(const VarContainer &var)
        {
            (*_val)(var);
        }

        double operator()() const
        {
            return (*_val)();
        }

    private:
        void construct(bool allocate = true)
        {
            this->_name = "pfunct_Mn3V";
            this->_var_name.clear();
            this->_var_name.push_back("n3");
            this->_var_description.clear();
            this->_var_description.push_back("concentration");
            
            if(!allocate) return;
            
            _val = new pfunct_Mn3V_f<VarContainer>();
        }

    };


}
#endif
