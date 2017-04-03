
#ifndef Hexahedron_HH
#define Hexahedron_HH

#include "../../pfunction/PSimpleBase.hh"
#include "../../pfunction/PSimpleFunction.hh"
#include "../../pfunction/PFuncBase.hh"
#include "./Interpolator.hh"
#include "../Coordinate.hh"

namespace PRISMS
{
    class Hexahedron_f : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            // var[0]: coordinate to be evaluated
            // var[1]: nodal coordinate
            // var[2]: element dimension
            // var[3]: +/- 1 depending on which 'corner' of quad
            //
            // f = (1.0 - e0)*(1.0 - e1)*(1.0 - e2)
            // e = var[3]*(var[0] - var[1])/var[2]
            
            return (1.0 - var[3][0]*(var[0][0] - var[1][0])/var[2][0])*
                   (1.0 - var[3][1]*(var[0][1] - var[1][1])/var[2][1])*
                   (1.0 - var[3][2]*(var[0][2] - var[1][2])/var[2][2]);
        }

    public:

        Hexahedron_f()
        {
            this->_name = "Hexahedron_f";
        }

        Hexahedron_f* clone() const
        {
            return new Hexahedron_f(*this);
        }
    };
    
    class Hexahedron_grad_0 : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            return -var[3][0]*(1.0 - var[3][1]*(var[0][1] - var[1][1])/var[2][1])*
                              (1.0 - var[3][2]*(var[0][2] - var[1][2])/var[2][2])/var[2][0];
        }

    public:

        Hexahedron_grad_0()
        {
            this->_name = "Hexahedron_grad_0";
        }

        Hexahedron_grad_0* clone() const
        {
            return new Hexahedron_grad_0(*this);
        }
    };
    
    class Hexahedron_grad_1 : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            return -var[3][1]*(1.0 - var[3][0]*(var[0][0] - var[1][0])/var[2][0])*
                              (1.0 - var[3][2]*(var[0][2] - var[1][2])/var[2][2])/var[2][1];
        
        }

    public:

        Hexahedron_grad_1()
        {
            this->_name = "Hexahedron_grad_1";
        }

        Hexahedron_grad_1* clone() const
        {
            return new Hexahedron_grad_1(*this);
        }
    };
    
    class Hexahedron_grad_2 : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            return -var[3][2]*(1.0 - var[3][0]*(var[0][0] - var[1][0])/var[2][0])*
                              (1.0 - var[3][1]*(var[0][1] - var[1][1])/var[2][1])/var[2][2];
        
        }

    public:

        Hexahedron_grad_2()
        {
            this->_name = "Hexahedron_grad_2";
        }

        Hexahedron_grad_2* clone() const
        {
            return new Hexahedron_grad_2(*this);
        }
    };
    
    class Hexahedron_hess_0_0 : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            return 0.0;
        }

    public:

        Hexahedron_hess_0_0()
        {
            this->_name = "Hexahedron_hess_0_0";
        }

        Hexahedron_hess_0_0* clone() const
        {
            return new Hexahedron_hess_0_0(*this);
        }
    };
    
    class Hexahedron_hess_0_1 : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            return var[3][0]*var[3][1]/var[2][0]/var[2][1];
        }

    public:

        Hexahedron_hess_0_1()
        {
            this->_name = "Hexahedron_hess_0_1";
        }

        Hexahedron_hess_0_1* clone() const
        {
            return new Hexahedron_hess_0_1(*this);
        }
    };
    
    class Hexahedron_hess_0_2 : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            return var[3][0]*var[3][2]/var[2][0]/var[2][2];
        }

    public:

        Hexahedron_hess_0_2()
        {
            this->_name = "Hexahedron_hess_0_2";
        }

        Hexahedron_hess_0_2* clone() const
        {
            return new Hexahedron_hess_0_2(*this);
        }
    };
    
    class Hexahedron_hess_1_0 : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            return var[3][1]*var[3][0]/var[2][1]/var[2][0];
        }

    public:

        Hexahedron_hess_1_0()
        {
            this->_name = "Hexahedron_hess_1_0";
        }

        Hexahedron_hess_1_0* clone() const
        {
            return new Hexahedron_hess_1_0(*this);
        }
    };
    
    class Hexahedron_hess_1_1 : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            return 0.0;
        }

    public:

        Hexahedron_hess_1_1()
        {
            this->_name = "Hexahedron_hess_1_1";
        }

        Hexahedron_hess_1_1* clone() const
        {
            return new Hexahedron_hess_1_1(*this);
        }
    };
    
    class Hexahedron_hess_1_2 : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            return var[3][1]*var[3][2]/var[2][1]/var[2][2];
        }

    public:

        Hexahedron_hess_1_2()
        {
            this->_name = "Hexahedron_hess_1_2";
        }

        Hexahedron_hess_1_2* clone() const
        {
            return new Hexahedron_hess_1_2(*this);
        }
    };
    
    class Hexahedron_hess_2_0 : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            return var[3][2]*var[3][0]/var[2][2]/var[2][0];
        }

    public:

        Hexahedron_hess_2_0()
        {
            this->_name = "Hexahedron_hess_2_0";
        }

        Hexahedron_hess_2_0* clone() const
        {
            return new Hexahedron_hess_2_0(*this);
        }
    };
    
    class Hexahedron_hess_2_1 : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            return var[3][2]*var[3][1]/var[2][2]/var[2][1];
        }

    public:

        Hexahedron_hess_2_1()
        {
            this->_name = "Hexahedron_hess_2_1";
        }

        Hexahedron_hess_2_1* clone() const
        {
            return new Hexahedron_hess_2_1(*this);
        }
    };
    
    class Hexahedron_hess_2_2 : public PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>
    {
        double eval( const std::vector<PRISMS::Coordinate<3> > &var) const
        {
            return 0.0;
        }

    public:

        Hexahedron_hess_2_2()
        {
            this->_name = "Hexahedron_hess_2_2";
        }

        Hexahedron_hess_2_2* clone() const
        {
            return new Hexahedron_hess_2_2(*this);
        }
    };
    
    
    class Hexahedron : public PFuncBase<std::vector<PRISMS::Coordinate<3> >, double>
    {
        PSimpleBase<std::vector<PRISMS::Coordinate<3> >, double>* _val;
        PSimpleBase<std::vector<PRISMS::Coordinate<3> >, double>** _grad_val;
        PSimpleBase<std::vector<PRISMS::Coordinate<3> >, double>*** _hess_val;
    
    public:    
        Hexahedron()
        {
            construct();
        }
        
        Hexahedron( const Hexahedron &RHS)
        {
            construct();
        }
        
        ~Hexahedron()
        {
            delete _val;

            delete _grad_val[0];
            delete _grad_val[1];
            delete _grad_val[2];
            delete [] _grad_val;

            delete _hess_val[0][0];
            delete _hess_val[0][1];
            delete _hess_val[0][2];
            delete _hess_val[1][0];
            delete _hess_val[1][1];
            delete _hess_val[1][2];
            delete _hess_val[2][0];
            delete _hess_val[2][1];
            delete _hess_val[2][2];
            delete [] _hess_val[0];
            delete [] _hess_val[1];
            delete [] _hess_val[2];
            delete [] _hess_val;
        }
        
        Hexahedron* clone() const
        {
            return new Hexahedron(*this);
        }

        PSimpleFunction< std::vector<PRISMS::Coordinate<3> >, double> simplefunction() const
        {
            return PSimpleFunction< std::vector<PRISMS::Coordinate<3> >, double>( *_val );
        }

        PSimpleFunction< std::vector<PRISMS::Coordinate<3> >, double> grad_simplefunction(size_type di) const
        {
            return PSimpleFunction< std::vector<PRISMS::Coordinate<3> >, double>( *_grad_val[di] );
        }

        PSimpleFunction< std::vector<PRISMS::Coordinate<3> >, double> hess_simplefunction(size_type di, size_type dj) const
        {
            return PSimpleFunction< std::vector<PRISMS::Coordinate<3> >, double>( *_hess_val[di][dj] );
        }

        double operator()(const std::vector<PRISMS::Coordinate<3> > &var)
        {
            return (*_val)(var);
        }

        double grad(const std::vector<PRISMS::Coordinate<3> > &var, size_type di)
        {
            return (*_grad_val[di])(var);
        }

        double hess(const std::vector<PRISMS::Coordinate<3> > &var, size_type di, size_type dj)
        {
            return (*_hess_val[di][dj])(var);
        }

        void eval(const std::vector<PRISMS::Coordinate<3> > &var)
        {
            (*_val)(var);
        }

        void eval_grad(const std::vector<PRISMS::Coordinate<3> > &var)
        {
            (*_grad_val[0])(var);
            (*_grad_val[1])(var);
        }

        void eval_hess(const std::vector<PRISMS::Coordinate<3> > &var)
        {
            (*_hess_val[0][0])(var);
            (*_hess_val[0][1])(var);
            (*_hess_val[1][0])(var);
            (*_hess_val[1][1])(var);
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
        void construct()
        {
            this->_name = "Hexahedron";
            this->_var_name.clear();
            this->_var_name.push_back("r");
            this->_var_name.push_back("n");
            this->_var_name.push_back("h");
            this->_var_name.push_back("s");
            this->_var_description.clear();
            this->_var_description.push_back("Coordinate to be evaluated (Cartesian)");
            this->_var_description.push_back("Coordinate of node");
            this->_var_description.push_back("Coordinate containing element dimensions");
            this->_var_description.push_back("Coordinate containing +/- 1.0, depending on which corner of quad element");
            
            _val = new Hexahedron_f();
            
            _grad_val = new PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>*[3];
            _grad_val[0] = new Hexahedron_grad_0();
            _grad_val[1] = new Hexahedron_grad_1();
            _grad_val[2] = new Hexahedron_grad_2();
            
            _hess_val = new PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>**[3];
            _hess_val[0] = new PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>*[3];
            _hess_val[1] = new PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>*[3];
            _hess_val[2] = new PSimpleBase< std::vector<PRISMS::Coordinate<3> >, double>*[3];
            _hess_val[0][0] = new Hexahedron_hess_0_0();
            _hess_val[0][1] = new Hexahedron_hess_0_1();
            _hess_val[0][2] = new Hexahedron_hess_0_2();
            _hess_val[1][0] = new Hexahedron_hess_1_0();
            _hess_val[1][1] = new Hexahedron_hess_1_1();
            _hess_val[1][2] = new Hexahedron_hess_1_2();
            _hess_val[2][0] = new Hexahedron_hess_2_0();
            _hess_val[2][1] = new Hexahedron_hess_2_1();
            _hess_val[2][2] = new Hexahedron_hess_2_2();
        }
    };
    
    
    /// A base class for interpolating functions
    ///
    template <class Coordinate>
    class HexahedronValues : public Interpolator<Coordinate, 3>
    {
        //_var[0]: Coordinate _r;  // coordinate to evaluate field at
        //_var[1]: Coordinate _n;  // coordinate of node
        //_var[2]: Coordinate _h;  // quad dimensions
        //_var[3]: Coordinate _s;  // +/- 1, depending on orientation of basis function
        std::vector< PRISMS::Coordinate<3> > _var;
        
    public:
        
        typedef typename Interpolator<Coordinate, 3>::size_type size_type;
    
        // node_index: index of node in mesh
        // node_index: index of element in mesh
        // bfunc: PFuncBase<std::vector<Coordinate>, double>* 
        // node_coord: Coordinate of node
        // dim: Coordinate containing x and y dimension of element
        // element_node_index: 0 == bottom left, proceed counter clockwise to 3 == top left of element
        
        HexahedronValues( unsigned long int node_index, 
                    unsigned long int element_index, 
                    PFuncBase<std::vector<PRISMS::Coordinate<3> >, double>* bfunc, 
                    const PRISMS::Coordinate<3> &node_coord, 
                    const PRISMS::Coordinate<3> &dim, 
                    int element_node_index)
        : Interpolator<Coordinate, 3>(node_index, element_index, bfunc)
        {
            _var.resize(4);
            
            _var[1][0] = node_coord[0];
            _var[1][1] = node_coord[1];
            _var[1][2] = node_coord[2];
                
            _var[2][0] = dim[0];
            _var[2][1] = dim[1];
            _var[2][2] = dim[2];
            
            if( element_node_index == 0)
            {
                _var[3][0] = 1.0;
                _var[3][1] = 1.0;
                _var[3][2] = 1.0;
            }
            else if( element_node_index == 1)
            {
                _var[3][0] = -1.0;
                _var[3][1] = 1.0;
                _var[3][2] = 1.0;
            }
            else if( element_node_index == 2)
            {
                _var[3][0] = -1.0;
                _var[3][1] = -1.0;
                _var[3][2] = 1.0;
            }
            else if( element_node_index == 3)
            {
                _var[3][0] = 1.0;
                _var[3][1] = -1.0;
                _var[3][2] = 1.0;
            }
            else if( element_node_index == 4)
            {
                _var[3][0] = 1.0;
                _var[3][1] = 1.0;
                _var[3][2] = -1.0;
            }
            else if( element_node_index == 5)
            {
                _var[3][0] = -1.0;
                _var[3][1] = 1.0;
                _var[3][2] = -1.0;
            }
            else if( element_node_index == 6)
            {
                _var[3][0] = -1.0;
                _var[3][1] = -1.0;
                _var[3][2] = -1.0;
            }
            else if( element_node_index == 7)
            {
                _var[3][0] = 1.0;
                _var[3][1] = -1.0;
                _var[3][2] = -1.0;
            }
        }
        
        PRISMS::Coordinate<3> min() const
        {
            PRISMS::Coordinate<3> coord = _var[1];
            
            if( _var[3][0] == -1.0)
                coord[0] -= _var[2][0];
            if( _var[3][1] == -1.0)
                coord[1] -= _var[2][1];
            if( _var[3][2] == -1.0)
                coord[2] -= _var[2][2];
            
            return coord;
        }
        
        PRISMS::Coordinate<3> max() const
        {
            PRISMS::Coordinate<3> coord = _var[1];
            
            if( _var[3][0] == 1.0)
                coord[0] += _var[2][0];
            if( _var[3][1] == 1.0)
                coord[1] += _var[2][1];
            if( _var[3][2] == 1.0)
                coord[2] += _var[2][2];
            
            return coord;
        }
    
        bool is_in_range(const Coordinate &coord)
        {
            _var[0][0] = coord[0];
            _var[0][1] = coord[1];
            _var[0][2] = coord[2];
            double e;
            
            for( int i=0; i<3; i++)
            {
                e = _var[3][i]*(_var[0][i] - _var[1][i])/_var[2][i];
                if( e < 0.0 || e >= 1.0)
                    return false;
                
                //if( e == 0.0 && std::signbit(e))
                //    return false;
            }
            
            
            
            //std::cout << "e: " ;
            //for( int i=0; i<2; i++)
            //{
            //    e = _var[3][i]*(_var[0][i] - _var[1][i])/_var[2][i];
            //    std::cout << e << " ";
            //}
            //std::cout << std::endl;
            
            return true;
        }
        
        // for the following,
        //   you are expected to KNOW that the coord is_in_range!!!
        
        double operator()(const Coordinate &coord)
        {
            _var[0][0] = coord[0];
            _var[0][1] = coord[1];
            _var[0][2] = coord[2];
            return (*this->_bfunc)(_var);
        }
        
        double grad(const Coordinate &coord, size_type di)
        {
            _var[0][0] = coord[0];
            _var[0][1] = coord[1];
            _var[0][2] = coord[2];
            return (*this->_bfunc).grad(_var,di);
        }
        
        double hess(const Coordinate &coord, size_type di, size_type dj)
        {
            _var[0][0] = coord[0];
            _var[0][1] = coord[1];
            _var[0][2] = coord[2];
            return (*this->_bfunc).hess(_var,di,dj);
        }

    };

}


#endif
