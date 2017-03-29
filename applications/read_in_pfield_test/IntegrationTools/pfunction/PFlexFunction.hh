
#ifndef PFlexFunction_HH
#define PFlexFunction_HH

#include<cstring>
#include<iostream>
#include<vector>
#include<cstdlib>

#include "./PFuncBase.hh"
#include "./PSimpleFunction.hh"

namespace PRISMS
{
    /// A class to create functions consisting of PSimpleFunctions
    ///  Used to create basis_functions from PSimpleFunctions for f, grad_f, and hess_f
    ///
    template<class VarContainer, class OutType>
    class PFlexFunction : public PFuncBase< VarContainer, OutType>
    {
        PSimpleFunction< VarContainer, OutType> _val;
        std::vector< PSimpleFunction< VarContainer, OutType> > _grad_val;
        std::vector< std::vector< PSimpleFunction< VarContainer, OutType> > > _hess_val;
        
    public:
        
        typedef typename PFuncBase< VarContainer, OutType>::size_type size_type;
    
        PFlexFunction()
        {
            
        }
        
        PFlexFunction( const std::string &name,
                  const std::vector<std::string> &var_name,
                  const std::vector<std::string> &var_description,
                  const PSimpleFunction< VarContainer, OutType> &simplef,
                  const std::vector< PSimpleFunction< VarContainer, OutType> > &grad_simplef,
                  const std::vector< std::vector< PSimpleFunction< VarContainer, OutType> > > &hess_simplef) : 
                  PFuncBase< VarContainer, OutType>(name, var_name, var_description), 
                  _val(simplef), _grad_val(grad_simplef), _hess_val(hess_simplef)
        {
            check();
        }
        
        void clear()
        {
            this->_name = "";
            this->_var_name.clear();
            this->_var_description.clear();
            _val = PSimpleFunction< VarContainer, OutType>();
            _grad_val.clear();
            _hess_val.clear();
        }
        
        void set( const std::string &name,
                  const std::vector<std::string> &var_name,
                  const std::vector<std::string> &var_description,
                  const PSimpleFunction< VarContainer, OutType> &simplef,
                  const std::vector< PSimpleFunction< VarContainer, OutType> > &grad_simplef,
                  const std::vector< std::vector< PSimpleFunction< VarContainer, OutType> > > &hess_simplef)
        {
            this->_name = name;
            this->_var_name = var_name;
            this->_var_description = var_description;
            _val = simplef;
            _grad_val = grad_simplef;
            _hess_val = hess_simplef;
            
            check();
        }

        PFlexFunction(const PFlexFunction &RHS )
        {
            this->_name = RHS._name;
            this->_var_name = RHS._var_name;
            this->_var_description = RHS._var_description;
            
            _val = RHS._val;
            _grad_val = RHS._grad_val;
            _hess_val = RHS._hess_val;
        }

        PFlexFunction& operator=(const PFlexFunction &RHS )
        {
            this->_name = RHS._name;
            this->_var_name = RHS._var_name;
            this->_var_description = RHS._var_description;
            
            _val = RHS._val;
            _grad_val = RHS._grad_val;
            _hess_val = RHS._hess_val;
            
        }

        ~PFlexFunction()
        {
            
        }

        PFlexFunction<VarContainer, OutType>* clone() const
        {
            return new PFlexFunction<VarContainer, OutType>(*this);
        }

        PSimpleFunction< VarContainer, OutType> simplefunction() const
        {
            return  _val;
        }

        PSimpleFunction< VarContainer, OutType> grad_simplefunction(size_type di) const
        {
            return _grad_val[di];
        }

        PSimpleFunction< VarContainer, OutType> hess_simplefunction(size_type di, size_type dj) const
        {
            return _hess_val[di][dj];
        }

        OutType operator()(const VarContainer &var)
        {
            return _val(var);
        }

        OutType grad(const VarContainer &var, size_type di)
        {
            return _grad_val[di](var);
        }

        OutType hess(const VarContainer &var, size_type di, size_type dj)
        {
            return _hess_val[di][dj](var);
        }
        
        void eval(const VarContainer &var)
        {
            (*this)(var);
        }

        void eval_grad(const VarContainer &var)
        {
            for( size_type i=0; i<_grad_val.size(); i++)
                _grad_val[i](var);
        }

        void eval_hess(const VarContainer &var)
        {
            for( size_type i=0; i<_hess_val.size(); i++)
                for( size_type j=0; j<_hess_val[i].size(); j++)
                    _hess_val[i][j](var);
        }

        OutType operator()() const
        {
            return _val();
        }

        OutType grad(size_type di) const
        {
            return _grad_val[di]();
        }

        OutType hess(size_type di, size_type dj) const
        {
            return _hess_val[di][dj]();
        }

        private:
        
        void check()
        {
            size_type n = this->_var_name.size();
            if( this->_var_description.size() != n)
            {
                std::cerr << "Error in PFlexFunction. _var_name.size() != _var_description.size()." << std::endl;
                exit(1);
            }
            if( _grad_val.size() != n)
            {
                std::cerr << "Error in PFlexFunction. _var_name.size() != _grad_val.size()." << std::endl;
                exit(1);
            }
            if( _hess_val.size() != n)
            {
                std::cerr << "Error in PFlexFunction. _var_name.size() != _hess_val.size()." << std::endl;
                exit(1);
            }
            for( size_type i=0; i<_hess_val.size(); i++)
            {
                if( _hess_val[i].size() != n)
                {
                    std::cerr << "Error in PFlexFunction. _var_name.size() != _hess_val[" << i << "].size()." << std::endl;
                    exit(1);
                }
            }
        }

    };
    
}


#endif