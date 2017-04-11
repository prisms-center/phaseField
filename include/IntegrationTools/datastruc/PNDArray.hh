
#ifndef PNDArray_HH
#define PNDArray_HH

#include <vector>
#include <iostream>
#include <cstdlib>

namespace PRISMS
{
    
    /// A Tensor class, 
    ///    takes int valued IndexContainer and returns OutType value
    ///
    template< class OutType>
    class PNDArray
    {
        std::vector< int > _dim;                        // dimension along each compenent
        std::vector< int > _unroll;                     // used to translate from tensor indices to linear index
        std::vector< OutType > _val;                    // unrolled list of coefficients (first index is outer loop)
        int _order;                                     // _dim.size()
        int _volume;                                    // _coeff_tensor.size() = product_i(_dim[i]) 
        
        public:
        
        PNDArray() : _order(0), _volume(0) {}
        PNDArray(const std::vector<int> &dim) 
        {
            resize(dim);
        }
        
        PNDArray(const std::vector<int> &dim, const std::vector<OutType> &value) 
        {
            resize(dim);
            if( _volume != value.size())
            {
                std::cerr << "Error in PNDArray(const std::vector<int> &dim, const std::vector<OutType> &value)." << std::endl;
                std::cerr << "  value.size() does not match volume based on dim." << std::endl;
                exit(1);
            }
            _val = value;
        }
        
        int order() const
        {
            return _order;
        }
        
        int volume() const
        {
            return _volume;
        }
        
        void resize(const std::vector<int> &dim)
        {
            _dim = dim;
            _order = _dim.size();
            _volume = calc_volume(_dim);
            _val.resize(_volume);
            generate_unroll();
        }
        
        void reshape(const std::vector<int> &dim)
        {
            if( calc_volume(dim) == _volume)
            {
                _dim = dim;
                _order = _dim.size();
                generate_unroll();
            }
            else
            {
                std::cerr << "Error in PNDArray::reshape. Volume is not equivalent." << std::endl;
                exit(1);
            }
        }
        
        void clear()
        {
            _val.clear();
            _dim.clear();
            _unroll.clear();
            _order = 0;
            _volume = 0;
            
        }
        
        const std::vector<int>& dim() const
        {
            return _dim;
        }
        
        int dim(int i) const
        {
            return _dim[i];
        }
        
        OutType& operator()( int i)
        {
            return _val[i];
        }
        
        template<class IndexContainer>
        OutType& operator()( const IndexContainer &term)
        {
            return _val[linear_index(term)];
        }
        
        template<class IndexContainer>
        int linear_index( const IndexContainer &term) const
        {
            int lindex = 0;
            for( unsigned int i=0; i<_unroll.size(); i++)
                lindex += term[i]*_unroll[i];
            return lindex;
        }
        
        template <class IndexContainer> 
        void tensor_indices( int lindex, IndexContainer &term) const
        {   // assumes term.size() == order()  (the tensor order)
            //   not sure if this is how we want to do it, but it avoids assuming push_back() or resize()
            
            for( int i=0; i<_unroll.size(); i++)
            {
                term[i] = lindex/_unroll[i];
                lindex -= term[i]*_unroll[i];
            }
        }
        
        
        private:
        
        int calc_volume( const std::vector<int> &dim)
        {
            if( dim.size() == 0)
            {
                return 0;
            }
            
            int vol = 1;
            for( unsigned int i=0; i<dim.size(); i++)
            {
                vol *= dim[i];
            }
            return vol;
        }
        
        void generate_unroll()
        {
            _unroll.resize(_dim.size());
            _unroll[_dim.size()-1] = 1;
            for( int i=_dim.size()-2; i>=0; i--)
                _unroll[i] = _unroll[i+1]*_dim[i+1];
        }
        
    };
}


#endif
