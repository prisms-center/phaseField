
#ifndef Mesh_HH
#define Mesh_HH

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cstdlib>

#include "../datastruc/Bin.hh"
#include "../pfunction/PFuncBase.hh"
#include "./interpolation/Interpolator.hh"
#include "./interpolation/Quad.hh"
#include "./interpolation/Hexahedron.hh"

namespace PRISMS
{

    inline void construct_basis_function( PFuncBase<std::vector<PRISMS::Coordinate<2> >, double>* &bfunc, const std::string &name)
    {
        if( name == "Quad")
        {
            bfunc = new Quad();
        }
        else
        {
            std::cout << "Error in construct_basis_function (2D): unknown name: " << name << std::endl;
            exit(1);
        }
    }

    inline void construct_basis_function( PFuncBase<std::vector<PRISMS::Coordinate<3> >, double>* &bfunc, const std::string &name)
    {
        if( name == "Hexahedron")
        {
            bfunc = new Hexahedron();
        }
        else
        {
            std::cout << "Error in construct_basis_function (3D): unknown name: " << name << std::endl;
            exit(1);
        }
    }

    template< class Coordinate>
    void construct_interpolating_functions( std::vector<Interpolator<Coordinate, 2>* > &interp,
                                        const std::string &name,
                                        unsigned long int cell,
                                        PFuncBase<std::vector<PRISMS::Coordinate<2> >, double>* bfunc_ptr,
                                        const std::vector<unsigned long int> &cell_node,
                                        const std::vector<PRISMS::Coordinate<2> > &node)
    {
        if( name == "Quad")
        {
            Interpolator<Coordinate, 2>* interp_ptr;

            //std::cout << "cell nodes: " << cell_node[0] << " " << cell_node[2] << std::endl;


            PRISMS::Coordinate<2> dim;
            dim[0] = node[ cell_node[2]][0] - node[ cell_node[0]][0];
            dim[1] = node[ cell_node[2]][1] - node[ cell_node[0]][1];

            // QuadValues(const Coordinate &node, const Coordinate &dim, int node_index)
            for( int j=0; j<4; j++)
            {
                interp.push_back(interp_ptr);
                interp.back() = new PRISMS::QuadValues<Coordinate>(cell_node[j], cell, bfunc_ptr, node[ cell_node[j]], dim, j);
            }
        }
        else
        {
            std::cout << "Error in construct_interpolating_function (2D): unknown name: " << name << std::endl;
            exit(1);
        }
    }

    template< class Coordinate>
    void construct_interpolating_functions( std::vector<Interpolator<Coordinate, 3>* > &interp,
                                        const std::string &name,
                                        unsigned long int cell,
                                        PFuncBase<std::vector<PRISMS::Coordinate<3> >, double>* bfunc_ptr,
                                        const std::vector<unsigned long int> &cell_node,
                                        const std::vector<PRISMS::Coordinate<3> > &node)
    {
        if( name == "Hexahedron")
        {
            Interpolator<Coordinate, 3>* interp_ptr;

            PRISMS::Coordinate<3> dim;
            dim[0] = node[ cell_node[6]][0] - node[ cell_node[0]][0];
            dim[1] = node[ cell_node[6]][1] - node[ cell_node[0]][1];
            dim[2] = node[ cell_node[6]][2] - node[ cell_node[0]][2];

            // QuadValues(const Coordinate &node, const Coordinate &dim, int node_index)
            for( int j=0; j<8; j++)
            {
                interp.push_back(interp_ptr);
                interp.back() = new PRISMS::HexahedronValues<Coordinate>(cell_node[j], cell, bfunc_ptr, node[ cell_node[j]], dim, j);
            }
        }
        else
        {
            std::cout << "Error in construct_interpolating_function (3D): unknown name: " << name << std::endl;
            exit(1);
        }
    }



    /// A template class for a finite element mesh
    ///   Needs: Coordinate::operator[]() for use in Bin
    ///
    template<class Coordinate, int DIM>
    class Mesh
    {
        // min and max coordinate of cuboid surrounding the body
        PRISMS::Coordinate<DIM> _min;
        PRISMS::Coordinate<DIM> _max;

        /// Vector of nodal coordinates
        ///    nodal values live in 'Field' class
        ///
        std::vector<PRISMS::Coordinate<DIM> > _node;

        /// array containing interpolating functions:
        ///    owns the interpolating functions
        ///    interpolating functions contain basis function / element info,
        ///    these point to _bfunc pfunctions which are used to evaluate
        ///
        std::vector<Interpolator<Coordinate, DIM>* > _interp;

        /// array containing PFunctions evaluated by interpolating functions
        ///    owns the pfunctions, which are pointed to by the interpolating functions
        ///    !!! do not modify after initial construction or pointers will be messed up !!!
        ///
        std::vector< PFuncBase<std::vector<PRISMS::Coordinate<DIM> >, double>* > _bfunc;

        /// bin of interpolating functions (this might be updated to be either Element or Spline Bins)
        ///
        Bin<Interpolator<Coordinate, DIM>*, Coordinate > _bin;

    public:

        // still need a constructor
        Mesh(){};

        ~Mesh()
        {
            for( unsigned int i=0; i<_interp.size(); i++)
            {
                delete _interp[i];
            }

            for( unsigned int i=0; i<_bfunc.size(); i++)
            {
                delete _bfunc[i];
            }
        };

        // reads vtk file through 'CELL_TYPES' and then returns
        void read_vtk(std::ifstream &infile)
        {
            bool mesh_as_points = true;
            std::vector<double> x_coord, y_coord, z_coord;

            std::istringstream ss;
            std::string line, str, type;

            unsigned long int uli_dummy;
            double d_dummy;

            unsigned long int Npoints, Ncells, Ncell_numbers;
            std::vector<unsigned long int> cell_node;

            PRISMS::Coordinate<DIM> _coord;

            while(!infile.eof())
            {
                std::getline( infile, line);
                //std::cout << "line: " << line << std::endl;

                if( line[0] == 'P')
                {
                    // read POINTS info:
                    // POINTS # type
                    // x y z
                    // x y z
                    // ...
                    if( line.size() > 5 && line.substr(0,6) == "POINTS")
                    {
                        // read header line
                        //std::cout << line << "\n";
                        ss.clear();
                        ss.str(line);
                        ss >> str >> Npoints >> type;

                        // read points
                        std::vector< std::vector<double> > value(DIM);
                        std::vector< std::vector<int> > hist(DIM);

                        std::cout << "Read POINTS: " << Npoints << std::endl;
                        _node.reserve(Npoints);
                        std::cout << "  reserve OK" << std::endl;
                        for( unsigned int i=0; i<Npoints; i++)
                        {
                            if( DIM == 2)
                            {
                                infile >> _coord[0] >> _coord[1] >> d_dummy;
                                //std::cout << _coord[0] << " " << _coord[1] << " " << d_dummy << std::endl;
                            }
                            else if( DIM == 3)
                            {
                                infile >> _coord[0] >> _coord[1] >> _coord[2];
                                //std::cout << _coord[0] << " " << _coord[1] << " " << _coord[3] << std::endl;
                            }

                            for( int j=0; j<DIM; j++)
                            add_once( value[j], hist[j], _coord[j]);
                            _node.push_back(_coord);
                        }
                        std::cout << "  done" << std::endl;

                        // create bins
                        std::vector<double> min;
                        std::vector<int> N;
                        std::vector<double> incr;

                        std::cout << "Determine Body size" << std::endl;
                        for( int j=0; j<DIM; j++)
                        {
                            for( unsigned int i=1; i<hist[j].size(); i++)
                            {
                                /*if( hist[j][i] != hist[j][i-1])
                                {
                                std::cout << "Error reading 2D vtk file." << std::endl;
                                std::cout << "  Not rectangular." << std::endl;
                                std::cout << "  Dimension: " << j << std::endl;
                                std::cout << value[j][i-1] << ": " << hist[j][i-1] << "  " << value[j][i] << ": " << hist[j][i] << std::endl;
                                exit(1);
                            }
                            */
                        }

                        std::sort( value[j].begin(), value[j].end());
                        //std::cout << "j: " << j << " back(): " << value[j].back() << std::endl;
                        min.push_back( value[j][0]);
                        N.push_back( value[j].size());
                        incr.push_back( (value[j].back() - value[j][0])/(1.0*N.back()));

                        // get min and max surrounding coordinates
                        _min[j] = value[j][0];
                        _max[j] = value[j].back();

                        // for short term, expand bin to avoid edge issues
                        min[j] -= incr[j];
                        N[j] += 2;
                    }
                    std::cout << "  Min Coordinate: ";
                    for( int j=0; j<DIM; j++)
                    std::cout << _min[j] << " ";
                    std::cout << std::endl;
                    std::cout << "  Max Coordinate: ";
                    for( int j=0; j<DIM; j++)
                    std::cout << _max[j] << " ";
                    std::cout << std::endl;

                    std::cout << "  done" << std::endl;

                    std::cout << "Initialize Bin" << std::endl;
                    _bin = Bin<Interpolator<Coordinate, DIM>*, Coordinate>(min, incr, N);
                    std::cout << "  done" << std::endl;

                }


            }
            else if( line[0] == 'C')
            {

                if( line.size() > 4 && line.substr(0,5) == "CELLS")
                {
                    //std::cout << line << "\n";
                    ss.clear();
                    ss.str(line);

                    ss >> str >> Ncells >> Ncell_numbers;

                    PFuncBase<std::vector<PRISMS::Coordinate<DIM> >, double>* bfunc_ptr;
                    _bfunc.push_back( bfunc_ptr);

                    if( DIM == 2)
                    {
                        // add Quad basis function
                        _interp.reserve(Ncells*4);
                        construct_basis_function(_bfunc.back(), "Quad");
                    }
                    else if( DIM == 3)
                    {
                        // add Hexahedron basis function
                        _interp.reserve(Ncells*8);
                        construct_basis_function(_bfunc.back(), "Hexahedron");
                    }
                    bfunc_ptr = _bfunc.back();

                    std::cout << "Read CELLS: " << Ncells << std::endl;
                    for( unsigned int i=0; i<Ncells; i++)
                    {
                        infile >> uli_dummy;

                        cell_node.resize(uli_dummy);
                        for( unsigned int j=0; j<uli_dummy; j++)
                        {
                            infile >> cell_node[j];
                        }

                        //std::cout << cell_node[0] << " " << cell_node[1] << " " << cell_node[2] << " " << cell_node[3] << std::endl;

                        // create interpolator
                        if( DIM == 2)
                        {
                            construct_interpolating_functions(_interp, "Quad", i, bfunc_ptr, cell_node, _node);
                        }
                        else if(DIM == 3)
                        {
                            construct_interpolating_functions(_interp, "Hexahedron", i, bfunc_ptr, cell_node, _node);
                        }

                    }
                    std::cout << "  done" << std::endl;

                    // bin interpolators
                    std::cout << "Bin interpolating functions" << std::endl;

                    for( unsigned int i=0; i<_interp.size(); i++)
                    {
                        _bin.add_range(_interp[i], _interp[i]->min(), _interp[i]->max());
                    }
                    std::cout << "  done  max_bin_size: " << _bin.max_size() << std::endl;


                }
                else if( line.size() > 9 && line.substr(0,10) == "CELL_TYPES")
                {
                    //std::cout << line << "\n";
                    ss.clear();
                    ss.str(line);

                    //std::cout << "ss.str()" << ss.str() << std::endl;
                    ss >> str >> Ncells;

                    for( unsigned int i=0; i<Ncells; i++)
                    {
                        infile >> uli_dummy;

                        if( uli_dummy != 9 && uli_dummy != 12)
                        {
                            std::cout << "Error reading CELL_TYPES: CELL TYPE != 9 && != 12" << std::endl;
                            std::cout << "   CELL TYPE: " << uli_dummy << std::endl;
                            exit(1);
                        }
                    }

                    return;
                }

            }
            if( line[0] == 'X')
            {
                mesh_as_points = false;

                // read X_COORDINATES info:
                // X_COORDINATES # type
                // x
                // x
                // ...
                if( line.size() > 12 && line.substr(0,13) == "X_COORDINATES")
                {
                    // read header line
                    //std::cout << line << "\n";
                    ss.clear();
                    ss.str(line);
                    ss >> str >> Npoints >> type;

                    std::cout << "Read X_COORDINATES: " << Npoints << std::endl;
                    _node.reserve(Npoints);
                    std::cout << "  reserve OK" << std::endl;
                    for( unsigned int i=0; i<Npoints; i++)
                    {
                        double temp_coord;

                        infile >> temp_coord;

                        x_coord.push_back(temp_coord);

                    }
                }
            }
            if( line[0] == 'Y')
            {
                mesh_as_points = false;

                // read Y_COORDINATES info:
                // Y_COORDINATES # type
                // y
                // y
                // ...
                if( line.size() > 12 && line.substr(0,13) == "Y_COORDINATES")
                {


                    // read header line
                    //std::cout << line << "\n";
                    ss.clear();
                    ss.str(line);
                    ss >> str >> Npoints >> type;

                    // read points
                    std::vector< std::vector<double> > value(DIM);
                    std::vector< std::vector<int> > hist(DIM);

                    std::cout << "Read Y_COORDINATES: " << Npoints << std::endl;
                    _node.reserve(Npoints);
                    std::cout << "  reserve OK" << std::endl;
                    for( unsigned int i=0; i<Npoints; i++)
                    {
                        double temp_coord;

                        infile >> temp_coord;

                        y_coord.push_back(temp_coord);

                    }
                }
            }
            if( line[0] == 'Z')
            {
                mesh_as_points = false;

                // read Z_COORDINATES info:
                // Z_COORDINATES # type
                // z
                // z
                // ...
                if( line.size() > 12 && line.substr(0,13) == "Z_COORDINATES")
                {


                    // read header line
                    //std::cout << line << "\n";
                    ss.clear();
                    ss.str(line);
                    ss >> str >> Npoints >> type;

                    // read points

                    std::cout << "Read Z_COORDINATES: " << Npoints << std::endl;
                    _node.reserve(Npoints);
                    std::cout << "  reserve OK" << std::endl;
                    for( unsigned int i=0; i<Npoints; i++)
                    {
                        double temp_coord;

                        infile >> temp_coord;

                        z_coord.push_back(temp_coord);

                    }
                }
            }
        }

        if (!mesh_as_points){
            std::vector< std::vector<double> > value(DIM);
            std::vector< std::vector<int> > hist(DIM);

            for (unsigned int i=0; i<x_coord.size(); i++){
                for (unsigned int j=0; j<y_coord.size(); j++){
                    for (unsigned int k=0; k<z_coord.size(); k++){
                        _coord[0] = x_coord.at(i);
                        _coord[1] = y_coord.at(j);
                        if (DIM > 2){
                            _coord[2] = z_coord.at(k);
                        }

                        for( int m=0; m<DIM; m++)
                        add_once( value[m], hist[m], _coord[m]);

                        _node.push_back(_coord);
                    }
                }
            }


            std::cout << "  done" << std::endl;

            // create bins
            std::vector<double> min;
            std::vector<int> N;
            std::vector<double> incr;

            std::cout << "Determine Body size" << std::endl;
            for( int j=0; j<DIM; j++)
            {
                for( unsigned int i=1; i<hist[j].size(); i++)
                {
                    /*if( hist[j][i] != hist[j][i-1])
                    {
                    std::cout << "Error reading 2D vtk file." << std::endl;
                    std::cout << "  Not rectangular." << std::endl;
                    std::cout << "  Dimension: " << j << std::endl;
                    std::cout << value[j][i-1] << ": " << hist[j][i-1] << "  " << value[j][i] << ": " << hist[j][i] << std::endl;
                    exit(1);
                }
                */
            }

            std::sort( value[j].begin(), value[j].end());
            //std::cout << "j: " << j << " back(): " << value[j].back() << std::endl;
            min.push_back( value[j][0]);
            N.push_back( value[j].size());
            incr.push_back( (value[j].back() - value[j][0])/(1.0*N.back()));

            // get min and max surrounding coordinates
            _min[j] = value[j][0];
            _max[j] = value[j].back();

            // for short term, expand bin to avoid edge issues
            min[j] -= incr[j];
            N[j] += 2;
        }
        std::cout << "  Min Coordinate: ";
        for( int j=0; j<DIM; j++)
        std::cout << _min[j] << " ";
        std::cout << std::endl;
        std::cout << "  Max Coordinate: ";
        for( int j=0; j<DIM; j++)
        std::cout << _max[j] << " ";
        std::cout << std::endl;

        std::cout << "  done" << std::endl;

        std::cout << "Initialize Bin" << std::endl;
        _bin = Bin<Interpolator<Coordinate, DIM>*, Coordinate>(min, incr, N);
        std::cout << "  done" << std::endl;

        // Now add the cell data
        unsigned int Ncells = (x_coord.size()-1) * (y_coord.size()-1);
        if (DIM > 2){
            Ncells *= (z_coord.size()-1);
        }

        PFuncBase<std::vector<PRISMS::Coordinate<DIM> >, double>* bfunc_ptr;
        _bfunc.push_back( bfunc_ptr);

        if( DIM == 2)
        {
            // add Quad basis function
            _interp.reserve(Ncells*4);
            construct_basis_function(_bfunc.back(), "Quad");
        }
        else if( DIM == 3)
        {
            // add Hexahedron basis function
            _interp.reserve(Ncells*8);
            construct_basis_function(_bfunc.back(), "Hexahedron");
        }
        bfunc_ptr = _bfunc.back();

        std::cout << "Read CELLS: " << Ncells << std::endl;

        unsigned int uli_dummy;
        if (DIM > 2){
            uli_dummy = 8;
        }
        else {
            uli_dummy = 4;
        }
        for( unsigned int i=0; i<Ncells; i++)
        {
            cell_node.resize(uli_dummy);
            for( unsigned int j=0; j<uli_dummy; j++)
            {
                cell_node[j] = i*uli_dummy+j;
            }

            if( DIM == 2)
            {
                double temp = cell_node[2];
                cell_node[2] = cell_node[3];
                cell_node[3] = temp;
            }

            //std::cout << cell_node[0] << " " << cell_node[1] << " " << cell_node[2] << " " << cell_node[3] << std::endl;

            // create interpolator
            if( DIM == 2)
            {
                construct_interpolating_functions(_interp, "Quad", i, bfunc_ptr, cell_node, _node);
            }
            else if(DIM == 3)
            {
                construct_interpolating_functions(_interp, "Hexahedron", i, bfunc_ptr, cell_node, _node);
            }

        }
        std::cout << "  done" << std::endl;

        // bin interpolators
        std::cout << "Bin interpolating functions" << std::endl;
        std::cout << "num nodes: " << _node.size() << std::endl;
        for( unsigned int i=0; i<_interp.size(); i++)
        {
            std::cout << "interp: " << _interp[i] << " " << _interp[i]->min() << " " << _interp[i]->max() << std::endl;
            _bin.add_range(_interp[i], _interp[i]->min(), _interp[i]->max());
        }
        std::cout << "  done  max_bin_size: " << _bin.max_size() << std::endl;

    }


}

        void min( Coordinate &coord)
        {
            for( int i=0; i<DIM; i++)
                coord[i] = _min[i];
        }

        void max( Coordinate &coord)
        {
            for( int i=0; i<DIM; i++)
                coord[i] = _max[i];
        }

        double min( int i)
        {
            return _min[i];
        }

        double max( int i)
        {
            return _max[i];
        }


        int max_bin_size()
        {
            return _bin.max_size();
        }

        // Set 'bfunc' to evaluated basis functions at 'coord',
        //     'node_index' to node indices for each basis function,
        //     and 's' is the length (number of basis functions)
        //  - 'bfunc' and 'node_index' are not resized, they must be big enough
        //
        void basis_functions(const Coordinate &coord, std::vector<double> &bfunc, std::vector<unsigned long int> &node_index, int &s)
        {
            std::vector<Interpolator<Coordinate,DIM>* > &bin = _bin.contents(coord);
            s = bin.size();

            int i=0;
            unsigned long int element;

            for( i=0; i<s; i++)
            {
                if( (*bin[i]).is_in_range(coord))
                {
                    element = (*bin[i]).element();
                    for( i=0; i<s; i++)
                    {
                        if( (*bin[i]).element() == element)
                        {
                            bfunc[i] = (*bin[i])(coord);
                        }
                        else
                        {
                            bfunc[i] = 0.0;
                        }
                        node_index[i] = (*bin[i]).node();
                        //std::cout << "i: " << i << "  bfunc: " << bfunc[i] << "  node: " << _node[ node_index[i]] << std::endl;
                    }
                    return;
                }
                //else
                //{
                //    bfunc[i] = 0.0;
                //    node_index[i] = (*bin[i]).node();
                //}
                //std::cout << "i: " << i << "  bfunc: " << bfunc[i] << "  node: " << _node[ node_index[i]] << std::endl;
            }

        };

        // Set 'bfunc' to evaluated grad basis functions at coord, and 's' is the length
        void grad_basis_functions(const Coordinate &coord, int di, std::vector<double> &bfunc, std::vector<unsigned long int> &node_index, int &s)
        {
            //std::cout << "begin Mesh::grad_basis_functions()" << std::endl;
            std::vector<Interpolator<Coordinate,DIM>* > &bin = _bin.contents(coord);
            s = bin.size();

            int i=0;
            unsigned long int element;

            for( i=0; i<s; i++)
            {
                if( (*bin[i]).is_in_range(coord))
                {
                    element = (*bin[i]).element();
                    for( i=0; i<s; i++)
                    {
                        if( (*bin[i]).element() == element)
                        {
                            bfunc[i] = (*bin[i]).grad(coord, di);
                        }
                        else
                        {
                            bfunc[i] = 0.0;
                        }
                        node_index[i] = (*bin[i]).node();
                        //std::cout << "i: " << i << "  bfunc: " << bfunc[i] << "  node: " << _node[ node_index[i]] << std::endl;
                    }
                    return;
                }
                //else
                //{
                //    bfunc[i] = 0.0;
                //    node_index[i] = (*bin[i]).node();
                //}
                //std::cout << "i: " << i << "  bfunc: " << bfunc[i] << "  node: " << _node[ node_index[i]] << std::endl;
            }
            //std::cout << "finish Mesh::grad_basis_functions()" << std::endl;

        }

        // Set 'bfunc' to evaluated hess basis functions at coord, and 's' is the length
        void hess_basis_functions(Coordinate coord, int di, int dj, std::vector<double> &bfunc, std::vector<unsigned long int> &node_index, int &s)
        {
            std::vector<Interpolator<Coordinate,DIM>* > &bin = _bin.contents(coord);
            s = bin.size();

            int i=0;
            unsigned long int element;

            for( i=0; i<s; i++)
            {
                if( (*bin[i]).is_in_range(coord))
                {
                    element = (*bin[i]).element();
                    for( i=0; i<s; i++)
                    {
                        if( (*bin[i]).element() == element)
                        {
                            bfunc[i] = (*bin[i]).hess(coord, di, dj);
                        }
                        else
                        {
                            bfunc[i] = 0.0;
                        }
                        node_index[i] = (*bin[i]).node();
                        //std::cout << "i: " << i << "  bfunc: " << bfunc[i] << "  node: " << _node[ node_index[i]] << std::endl;
                    }
                    return;
                }
                //else
                //{
                //    bfunc[i] = 0.0;
                //    node_index[i] = (*bin[i]).node();
                //}
                //std::cout << "i: " << i << "  bfunc: " << bfunc[i] << "  node: " << _node[ node_index[i]] << std::endl;
            }
        }

    private:

        void add_once( std::vector<double> &list, std::vector<int> &hist, double val)
        {
            //std::cout << "begin add_once()" << std::endl;

            for( unsigned int i=0; i<list.size(); i++)
            {
                if( list[i] == val)
                {
                    hist[i]++;
                    return;
                }
            }

            list.push_back(val);
            hist.push_back(1);

            //std::cout << "finish add_once()" << std::endl;

        }
    };

}


#endif
