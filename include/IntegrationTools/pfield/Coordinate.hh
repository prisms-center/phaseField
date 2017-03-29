
#ifndef Coordinate_HH
#define Coordinate_HH

namespace PRISMS
{
    
    /// A class for a coordinate, templated by dimension
    ///    This is a possible option anyplace 'Coordinate' class template is used
    ///    but it is not the only option. Any class that implements
    ///    'Coordinate::operator[]()' should work
    
    template<int DIM>
    class Coordinate
    {
        double _coord[DIM];
    public:
    
        
        int size() const
        {
            return DIM;
        }
        
        double& operator[](int i)
        {
            return _coord[i];
        }
        
        const double& operator[](int i) const
        {
            return _coord[i];
        }
        
        template< int D> friend std::ostream &operator<<(std::ostream &outstream, const Coordinate<D> &coord);
        
    };
    
    template<int D> std::ostream &operator<<(std::ostream &outstream, const Coordinate<D> &coord)
    {
        for( int i=0; i<coord.size(); i++)
        {
            outstream << coord[i];
            if( i < coord.size()-1)
                outstream << " ";
        }
        return outstream;
    }
}


#endif