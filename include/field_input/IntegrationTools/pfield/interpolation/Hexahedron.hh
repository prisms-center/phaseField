#ifndef Hexahedron_HH
#define Hexahedron_HH

#include <field_input/IntegrationTools/pfield/Coordinate.hh>
#include <field_input/IntegrationTools/pfield/interpolation/Interpolator.hh>
#include <field_input/IntegrationTools/pfunction/PFuncBase.hh>
#include <field_input/IntegrationTools/pfunction/PSimpleBase.hh>

namespace PRISMS
{
  class Hexahedron_f : public PSimpleBase<std::vector<PRISMS::Coordinate<3>>, double>
  {
    [[nodiscard]] double
    eval(const std::vector<PRISMS::Coordinate<3>> &var) const override
    {
      // var[0]: coordinate to be evaluated
      // var[1]: nodal coordinate
      // var[2]: element dimension
      // var[3]: +/- 1 depending on which 'corner' of quad
      //
      // f = (1.0 - e0)*(1.0 - e1)*(1.0 - e2)
      // e = var[3]*(var[0] - var[1])/var[2]

      return (1.0 - var[3][0] * (var[0][0] - var[1][0]) / var[2][0]) *
             (1.0 - var[3][1] * (var[0][1] - var[1][1]) / var[2][1]) *
             (1.0 - var[3][2] * (var[0][2] - var[1][2]) / var[2][2]);
    }

  public:
    Hexahedron_f()
    {
      this->_name = "Hexahedron_f";
    }
  };

  class Hexahedron : public PFuncBase<std::vector<PRISMS::Coordinate<3>>, double>
  {
    PSimpleBase<std::vector<PRISMS::Coordinate<3>>, double> *_val;

  public:
    Hexahedron()
      : PFuncBase<std::vector<PRISMS::Coordinate<3>>, double>()
    {
      _val = new Hexahedron_f();
    }

    ~Hexahedron()
    {
      delete _val;
    }

    double
    operator()(const std::vector<PRISMS::Coordinate<3>> &var) override
    {
      return (*_val)(var);
    }

    void
    eval(const std::vector<PRISMS::Coordinate<3>> &var) override
    {
      (*_val)(var);
    }

    double
    operator()() const override
    {
      return (*_val)();
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
    std::vector<PRISMS::Coordinate<3>> _var;

  public:
    using size_type = typename Interpolator<Coordinate, 3>::size_type;

    // node_index: index of node in mesh
    // node_index: index of element in mesh
    // bfunc: PFuncBase<std::vector<Coordinate>, double>*
    // node_coord: Coordinate of node
    // dim: Coordinate containing x and y dimension of element
    // element_node_index: 0 == bottom left, proceed counter clockwise to 3 == top left of
    // element

    HexahedronValues(unsigned long int                                      node_index,
                     unsigned long int                                      element_index,
                     PFuncBase<std::vector<PRISMS::Coordinate<3>>, double> *bfunc,
                     const PRISMS::Coordinate<3>                           &node_coord,
                     const PRISMS::Coordinate<3>                           &dim,
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

      if (element_node_index == 0)
        {
          _var[3][0] = 1.0;
          _var[3][1] = 1.0;
          _var[3][2] = 1.0;
        }
      else if (element_node_index == 1)
        {
          _var[3][0] = -1.0;
          _var[3][1] = 1.0;
          _var[3][2] = 1.0;
        }
      else if (element_node_index == 2)
        {
          _var[3][0] = -1.0;
          _var[3][1] = -1.0;
          _var[3][2] = 1.0;
        }
      else if (element_node_index == 3)
        {
          _var[3][0] = 1.0;
          _var[3][1] = -1.0;
          _var[3][2] = 1.0;
        }
      else if (element_node_index == 4)
        {
          _var[3][0] = 1.0;
          _var[3][1] = 1.0;
          _var[3][2] = -1.0;
        }
      else if (element_node_index == 5)
        {
          _var[3][0] = -1.0;
          _var[3][1] = 1.0;
          _var[3][2] = -1.0;
        }
      else if (element_node_index == 6)
        {
          _var[3][0] = -1.0;
          _var[3][1] = -1.0;
          _var[3][2] = -1.0;
        }
      else if (element_node_index == 7)
        {
          _var[3][0] = 1.0;
          _var[3][1] = -1.0;
          _var[3][2] = -1.0;
        }
    }

    [[nodiscard]] PRISMS::Coordinate<3>
    min() const override
    {
      PRISMS::Coordinate<3> coord = _var[1];

      if (_var[3][0] == -1.0)
        {
          coord[0] -= _var[2][0];
        }
      if (_var[3][1] == -1.0)
        {
          coord[1] -= _var[2][1];
        }
      if (_var[3][2] == -1.0)
        {
          coord[2] -= _var[2][2];
        }

      return coord;
    }

    [[nodiscard]] PRISMS::Coordinate<3>
    max() const override
    {
      PRISMS::Coordinate<3> coord = _var[1];

      if (_var[3][0] == 1.0)
        {
          coord[0] += _var[2][0];
        }
      if (_var[3][1] == 1.0)
        {
          coord[1] += _var[2][1];
        }
      if (_var[3][2] == 1.0)
        {
          coord[2] += _var[2][2];
        }

      return coord;
    }

    bool
    is_in_range(const Coordinate &coord) override
    {
      _var[0][0] = coord[0];
      _var[0][1] = coord[1];
      _var[0][2] = coord[2];
      double e;

      for (int i = 0; i < 3; i++)
        {
          e = _var[3][i] * (_var[0][i] - _var[1][i]) / _var[2][i];
          if (e < 0.0 || e >= 1.0)
            {
              return false;
            }
        }

      return true;
    }

    // for the following,
    //   you are expected to KNOW that the coord is_in_range!!!

    double
    operator()(const Coordinate &coord) override
    {
      _var[0][0] = coord[0];
      _var[0][1] = coord[1];
      _var[0][2] = coord[2];
      return (*this->_bfunc)(_var);
    }
  };

} // namespace PRISMS

#endif
