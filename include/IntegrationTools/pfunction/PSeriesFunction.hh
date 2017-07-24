
#ifndef PSeriesFunction_HH
#define PSeriesFunction_HH

#include<cstring>
#include<iostream>
#include<vector>
#include<limits>

#include "./PFuncBase.hh"
#include "./PBasisSet.hh"
#include "../datastruc/PNDArray.hh"

namespace PRISMS
{
    /// This class is for a series function that is the tensor product of indepedent bases
    ///
    ///   Example for order-3 tensor:
    ///
    ///    f(x,y,z) = sum_i sum_j sum_k T_ijk * phix_i(x) * phiy_j(y) * phiz_k(z)
    ///
    ///    summations run from i to _dim[i]
    ///    x, y, z are independent variables (dx/dy = 0, etc.)
    ///
    /// Template parameters:
    ///   - InType: input type for BasisSets (must be the same for all BasisSets)
    ///   - OutType: output type for SeriesFunction and BasisSets
    ///   - VarContainer: container for input
    ///   - IndexContainer: container for tensor indices
    ///
    /// Example template parameters are: <double, double, std::vector<double>, std::vector<int> >
    ///
    template< class InType, class OutType, class VarContainer, class IndexContainer>
    class PSeriesFunction : public PFuncBase<VarContainer, OutType>
    {
    protected:
        
        // coefficient tensor
        PNDArray<OutType> _coeff;
        
        OutType _identity_val;
        OutType _zero_val;
        
        // basis sets
        std::vector< PBasisSet<InType, OutType> > _basis_set;
        
        // evaluated values
        OutType _val;
        std::vector<OutType> _grad_val;
        std::vector< std::vector< OutType> > _hess_val;
        
    public:    
        PSeriesFunction(OutType zero, OutType identity)
        {
            _zero_val = zero;
            _identity_val = identity;
        }
        
        PSeriesFunction( OutType zero, OutType identity, const std::vector<PBasisSet<InType,OutType> > &basis_set) 
        {
            _zero_val = zero;
            _identity_val = identity;

            set(basis_set);
        }
        
        virtual ~PSeriesFunction()
        {
            clear();
        }
        
        void clear()
        {
            _basis_set.clear();
            
            _coeff.clear();
            
            _val = _zero_val;
            _grad_val.clear();
            _hess_val.clear();
            
        }
        
        /// generate the PSeriesFunction, by cloning 'basis_set' and resizing everything to match
        ///
        void set( const std::vector<PBasisSet<InType,OutType> > &basis_set)
        {
            //std::cout << "begin PSeriesFunction::set()" << std::endl;
            clear();
            
            std::vector<int> dim(basis_set.size());
            for( int i=0; i<basis_set.size(); i++)
                dim[i] = basis_set[i].size();
            _coeff.resize(dim);
            
            //std::cout << "  resize: " << basis_set.size() << std::endl;
            _basis_set.resize(basis_set.size());
            
            _grad_val.resize( _coeff.order() );
            _hess_val.resize( _coeff.order() );
            
            for( int i=0; i<_coeff.order(); i++)
            {
                //std::cout << "  i: " << i << std::endl;
                _hess_val[i].resize( _coeff.order() );
                //std::cout << "  copy basis set" << std::endl;
                _basis_set[i] = basis_set[i];
            }
            
        }
        
        PNDArray<OutType>& coeff()
        {
            return _coeff;
        }

        virtual PSeriesFunction<InType, OutType, VarContainer, IndexContainer> *clone() const
        {
            return new PSeriesFunction<InType, OutType, VarContainer, IndexContainer>(*this);
        }

        // ----------------------------------------------------------
        // Use these functions if you want to evaluate a single value
        virtual OutType operator()(const VarContainer &var)
        {
            // evaluate basis functions needed
            eval_basis(var);
            
            // evaluate basis function products, multiply by _coeff_tensor & sum
            return calc_val();
        }
        
        virtual OutType grad(const VarContainer &var, int di)
        {
            // evaluate basis functions needed
            for( int i=0; i<_coeff.order(); i++)
            {
                if( i == di)
                {   
                    _basis_set[i].eval_grad(var[i]);
                }
                else
                {
                    _basis_set[i](var[i]);
                }
            }
            
            // evaluate basis function products, multiply by _coeff_tensor & sum
            return calc_grad_val(di);
        }
        
        virtual OutType hess(const VarContainer &var, int di, int dj)
        {
            // evaluate basis functions needed
            for( int i=0; i<_coeff.order(); i++)
            {
                if( i == di || i == dj)
                {   
                    if( di == dj)
                    {
                        _basis_set[i].eval_hess(var[i]);
                    }
                    else
                    {
                        _basis_set[i].eval_grad(var[i]);
                    }
                }
                else
                {
                    _basis_set[i](var[i]);
                }
            }
            
            // evaluate basis function products, multiply by _coeff_tensor & sum
            return calc_hess_val(di,dj);
        }
        
        // ----------------------------------------------------------
        // Use these functions to evaluate several values, then use 'get' methods to access results
        virtual void eval(const VarContainer &var)
        {
            (*this)(var);
        }
        virtual void eval_grad(const VarContainer &var)
        {
            eval_basis(var);
            
            eval_basis_grad(var);
            
            for( int i=0; i<_coeff.order(); i++)
                (*this).calc_grad_val(i);
        }
        virtual void eval_hess(const VarContainer &var)
        {
            eval_basis(var);
            
            eval_basis_grad(var);
            
            eval_basis_hess(var);
            
            
            for( int i=0; i<_coeff.order(); i++)
            for( int j=0; j<_coeff.order(); j++)
                (*this).calc_hess_val(i,j);
        }

        virtual OutType operator()() const
        {
            return _val;
        }
        virtual OutType grad(int di) const
        {
            return _grad_val[di];
        }
        virtual OutType hess(int di, int dj) const
        {
            return _hess_val[di][dj];
        }


        // Functions for evaluating basis functions & their derivatives:

        // ----------------------------------------------------------
        // Use these functions if you want to evaluate a single value

        //   use basis index and term index for individual basis function
        virtual OutType basis(int bindex, int term, const VarContainer &var)
        {
            return _basis_set[bindex](term,var[bindex]);
        }
        virtual OutType basis_grad(int bindex, int term, const VarContainer &var)
        {
            return _basis_set[bindex].grad(term, var[bindex]);
        }
        virtual OutType basis_hess(int bindex, int term, const VarContainer &var)
        {
            return _basis_set[bindex].hess(term, var[bindex]);
        }

        //   or use tensor indices to evaluate basis function product
        virtual OutType basis(const IndexContainer &term, const VarContainer &var)
        {
            // evaluate basis functions needed
            OutType tmp = _identity_val;
            for( int i=0; i<_coeff.order(); i++)
                tmp *= _basis_set[i](term[i],var[i]);
                
            return tmp;
        }
        virtual OutType basis_grad(const IndexContainer &term, const VarContainer &var, int di)
        {
            // evaluate basis functions needed
            OutType tmp = _identity_val;
            for( int i=0; i<_coeff.order(); i++)
            {
                if( i == di)
                {   
                    tmp *= _basis_set[i].grad(term[i], var[i]);
                }
                else
                {
                    tmp *= _basis_set[i](term[i], var[i]);
                }
            }
            
            return tmp;
        }
        virtual OutType basis_hess(const IndexContainer &term, const VarContainer &var, int di, int dj)
        {
            // evaluate basis functions needed
            OutType tmp = _identity_val;
            for( int i=0; i<_coeff.order(); i++)
            {
                if( i == di || i == dj)
                {   
                    if( di == dj)
                    {
                        tmp *= _basis_set[i].hess(term[i], var[i]);
                    }
                    else
                    {
                        tmp *= _basis_set[i].grad(term[i], var[i]);
                    }
                }
                else
                {
                    tmp *= _basis_set[i](term[i], var[i]);
                }
            }
            
            return tmp;
        }

        // ----------------------------------------------------------
        // Use these functions to evaluate all basis functions,
        //   then use following methods to access results.

        virtual void eval_basis(const VarContainer &var)
        {
            // evaluate basis functions
            for( int i=0; i<_coeff.order(); i++)
                _basis_set[i].eval(var[i]);
            
        }
        virtual void eval_basis(const VarContainer &var, int i)
        {
            // evaluate basis functions
            _basis_set[i].eval(var[i]);
            
        }
        virtual void eval_basis_grad(const VarContainer &var)
        {
            // evaluate basis grad functions
            for( int i=0; i<_coeff.order(); i++)
                _basis_set[i].eval_grad(var[i]);
        }
        virtual void eval_basis_grad(const VarContainer &var, int i)
        {
            // evaluate basis grad functions
            _basis_set[i].eval_grad(var[i]);
        }
        virtual void eval_basis_hess(const VarContainer &var)
        {
            // evaluate basis hess functions
            for( int i=0; i<_coeff.order(); i++)
                _basis_set[i].eval_hess(var[i]);
        }
        virtual void eval_basis_hess(const VarContainer &var, int i)
        {
            // evaluate basis hess functions
            _basis_set[i].eval_hess(var[i]);
        }

        //   use basis index and term index for individual basis function
        virtual OutType basis(int bindex, int term) const
        {
            return _basis_set[bindex](term);
        }
        virtual OutType basis_grad(int bindex, int term) const
        {
            return _basis_set[bindex].grad(term);
        }
        virtual OutType basis_hess(int bindex, int term) const
        {
            return _basis_set[bindex].hess(term);
        }

        //   or use tensor indices to evaluate basis function product
        virtual OutType basis(const IndexContainer &term) const
        {
            // evaluate basis function product
            OutType tmp = _identity_val;
            for( int i=0; i<_coeff.order(); i++)
                tmp *= _basis_set[i](term[i]);
                
            return tmp;
        }
        virtual OutType basis_grad(const IndexContainer &term, int di) const
        {
            // evaluate basis function product
            OutType tmp = _identity_val;
            for( int i=0; i<_coeff.order(); i++)
            {
                if( i == di)
                {   
                    tmp *= _basis_set[i].grad(term[i]);
                }
                else
                {
                    tmp *= _basis_set[i](term[i]);
                }
            }
            
            return tmp;
        }
        virtual OutType basis_hess(const IndexContainer &term, int di, int dj) const
        {
            // evaluate basis function product
            OutType tmp = _identity_val;
            for( int i=0; i<_coeff.order(); i++)
            {
                if( i == di || i == dj)
                {   
                    if( di == dj)
                    {
                        tmp *= _basis_set[i].hess(term[i]);
                    }
                    else
                    {
                        tmp *= _basis_set[i].grad(term[i]);
                    }
                }
                else
                {
                    tmp *= _basis_set[i](term[i]);
                }
            }
            
            return tmp;
        }
        
        // ----------------------------------------------------------
        // Read and write routines: 
        //   Write routines only uses 'getters', they assume you have 
        //   already evaluated the necessary basis functions
        //
        // These are not included in PExtern
        
        void print_basis(std::ostream &sout)
        {
            std::vector<int> term = std::vector<int>(_coeff.order(),0);
            for( int i=0; i<_coeff.volume(); i++)
            {
                _coeff.tensor_indices(i,term);
                for( int j=0; j<term.size(); j++)
                    sout << term[j] << " ";
                sout << basis(term) << "\n";
            }
                
        }
        void print_basis_grad(std::ostream &sout, int di)
        {
            std::vector<int> term = std::vector<int>(_coeff.order(),0);
            for( int i=0; i<_coeff.volume(); i++)
            {
                _coeff.tensor_indices(i,term);
                for( int j=0; j<term.size(); j++)
                    sout << term[j] << " ";
                sout << basis_grad(term,di) << "\n";
            }
                
        }
        void print_basis_hess(std::ostream &sout, int di, int dj)
        {
            std::vector<int> term = std::vector<int>(_coeff.order(),0);
            for( int i=0; i<_coeff.volume(); i++)
            {
                _coeff.tensor_indices(i,term);
                for( int j=0; j<term.size(); j++)
                    sout << term[j] << " ";
                sout << basis_hess(term,di,dj) << "\n";
            }
                
        }
        void print_coeff(std::ostream &sout)
        {
            std::vector<int> term = std::vector<int>(_coeff.order(),0);
            for( int i=0; i<_coeff.volume(); i++)
            {
                _coeff.tensor_indices(i,term);
                for( int j=0; j<term.size(); j++)
                    sout << term[j] << " ";
                sout << _coeff(i) << "\n";
            }
        }
        void read_coeff(std::istream &sin)
        {
            // no error checking, 
            //   assumes format identical to print_coeff output
            std::vector<int> term = std::vector<int>(_coeff.order(),0);
            for( int i=0; i<_coeff.volume(); i++)
            {
                // ignores tensor indices
                for( int j=0; j<term.size(); j++)
                {    
                    sin >> term[j];
                }
                sin >> _coeff(i);
            }
        }
        
        
        ///
        private:
        
        // ----------------------------------------------------------
        // These assume you've evaluated the necessary _basis, _basis_grad, and _basis_hess functions
        //   and just need to take products & sum
        
        // evaluate basis function products, multiply by _coeff_tensor & sum
        OutType calc_val()
        {
            std::vector<int> tindex(_coeff.order());
            OutType tmp;
            _val = _zero_val;
            for( int i=0; i<_coeff.volume(); i++)
            {
                tmp = _coeff(i);
                _coeff.tensor_indices(i, tindex);
                
                for( int j=0; j<_coeff.order(); j++)
                {
                    tmp *= _basis_set[j](tindex[j]);
                }
                _val += tmp;
            }
            
            return _val;
        }
        
        // evaluate basis function products, multiply by _coeff_tensor & sum
        OutType calc_grad_val(int di)
        {
            std::vector<int> tindex(_coeff.order());
            OutType tmp;
            _grad_val[di] = _zero_val;
            for( int i=0; i<_coeff.volume(); i++)
            {
                tmp = _coeff(i);
                _coeff.tensor_indices(i, tindex);
                for( int j=0; j<_coeff.order(); j++)
                {
                    if( j == di)
                    {
                        tmp *= _basis_set[j].grad(tindex[j]);
                    }
                    else
                    {
                        tmp *= _basis_set[j](tindex[j]);
                    }
                }
                
                _grad_val[di] += tmp;
                
            }
            return _grad_val[di];
        }
        
        // evaluate basis function products, multiply by _coeff_tensor & sum
        OutType calc_hess_val(int di, int dj)
        {
            std::vector<int> tindex(_coeff.order());
            OutType tmp;
            _hess_val[di][dj] = _zero_val;
            for( int i=0; i<_coeff.volume(); i++)
            {
                tmp = _coeff(i);
                _coeff.tensor_indices(i, tindex);
                for( int j=0; j<_coeff.order(); j++)
                {
                    if( j == di || j == dj)
                    {
                        if( di == dj)
                        {
                            tmp *= _basis_set[j].hess(tindex[j]);
                        }
                        else
                        {
                            tmp *= _basis_set[j].grad(tindex[j]);
                        }
                    }
                    else
                    {
                        tmp *= _basis_set[j](tindex[j]);
                    }
                }
                _hess_val[di][dj] += tmp;
            }
            return _hess_val[di][dj];
        }
    };
    

}


#endif