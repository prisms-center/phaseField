
#ifndef Interpolator_HH
#define Interpolator_HH

#include "../../pfunction/PFuncBase.hh"
#include "../Coordinate.hh"

namespace PRISMS
{
    
    /// A base class for interpolating functions
    ///
    template <class Coordinate, int DIM>
    class Interpolator
    {
    protected:
    
        unsigned long int _node;        //index of nodal value or control point
        unsigned long int _element;     //index of element
        PFuncBase<std::vector<PRISMS::Coordinate<DIM> >, double>* _bfunc;          // basis function to evaluate
        
    public:
        
        typedef typename PFuncBase<std::vector<PRISMS::Coordinate<DIM> >, double>::size_type size_type;
        
        Interpolator( unsigned long int node, unsigned long int element, PFuncBase<std::vector<PRISMS::Coordinate<DIM> >, double>* bfunc):
        _node(node), _element(element), _bfunc(bfunc)
        {};

        virtual ~Interpolator() {};
        
        unsigned long int node()
        {
            return _node;
        }
        
        unsigned long int element()
        {
            return _element;
        }
        
        virtual PRISMS::Coordinate<DIM> min() const
        {
            undefined("void min(Coordinate &coord) const");
            return PRISMS::Coordinate<DIM>();
        }
        
        virtual PRISMS::Coordinate<DIM> max() const
        {
            undefined("void max(Coordinate &coord) const");
            return PRISMS::Coordinate<DIM>();
        }
        
        virtual bool is_in_range(const Coordinate &coord)
        {
            undefined("bool is_in_range(Coordinate coord) const");
            return false;
        }
        
        virtual double operator()(const Coordinate &coord)
        {
            undefined("double operator()(Coordinate coord)");
            return double();
        }
        
        virtual double grad(const Coordinate &coord, size_type di)
        {
            undefined("double grad()(Coordinate coord, size_type di)");
            return double();
        }
        
        virtual double hess(const Coordinate &coord, size_type di, size_type dj)
        {
            undefined("double hess()(Coordinate coord, size_type di, size_type dj)");
            return double();
        }
        
    private:
        
        void undefined(std::string fname) const
        {
            std::cout << "Error in Interpolator." << std::endl;
            std::cout << "   The member function '" << fname << "' has not been defined." << std::endl;
            exit(1);
        }

    };

}


#endif
