#ifndef Quad_HH
#define Quad_HH

#include <field_input/IntegrationTools/pfield/Coordinate.hh>
#include <field_input/IntegrationTools/pfield/interpolation/Interpolator.hh>
#include <field_input/IntegrationTools/pfunction/PFuncBase.hh>
#include <field_input/IntegrationTools/pfunction/PSimpleBase.hh>

namespace PRISMS
{
  class Quad_f : public PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double>
  {
    [[nodiscard]] double
    eval(const std::vector<PRISMS::Coordinate<2>> &var) const override
    {
      // var[0]: coordinate to be evaluated
      // var[1]: nodal coordinate
      // var[2]: element dimension
      // var[3]: +/- 1 depending on which 'corner' of quad
      //
      // f = (1.0 - e0)*(1.0 - e1)
      // e = var[3]*(var[0] - var[1])/var[2]

      return (1.0 - var[3][0] * (var[0][0] - var[1][0]) / var[2][0]) *
             (1.0 - var[3][1] * (var[0][1] - var[1][1]) / var[2][1]);
    }

  public:
    Quad_f()
    {
      this->_name = "Quad_f";
    }
  };

  class Quad : public PFuncBase<std::vector<PRISMS::Coordinate<2>>, double>
  {
    PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double> *_val;

  public:
    using size_type = PFuncBase<std::vector<PRISMS::Coordinate<2>>, double>::size_type;

    Quad()
      : PFuncBase<std::vector<PRISMS::Coordinate<2>>, double>()
    {
      _val = new Quad_f();
    }

    ~Quad()
    {
      delete _val;
    }

    double
    operator()(const std::vector<PRISMS::Coordinate<2>> &var) override
    {
      return (*_val)(var);
    }

    void
    eval(const std::vector<PRISMS::Coordinate<2>> &var) override
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
  class QuadValues : public Interpolator<Coordinate, 2>
  {
    //_var[0]: Coordinate _r;  // coordinate to evaluate field at
    //_var[1]: Coordinate _n;  // coordinate of node
    //_var[2]: Coordinate _h;  // quad dimensions
    //_var[3]: Coordinate _s;  // +/- 1, depending on orientation of basis function
    std::vector<PRISMS::Coordinate<2>> _var;

  public:
    using size_type = typename Interpolator<Coordinate, 2>::size_type;

    // node_index: index of node in mesh
    // node_index: index of element in mesh
    // bfunc: PFuncBase<std::vector<Coordinate>, double>*
    // node_coord: Coordinate of node
    // dim: Coordinate containing x and y dimension of element
    // element_node_index: 0 == bottom left, proceed counter clockwise to 3 == top left of
    // element

    QuadValues(unsigned long int                                      node_index,
               unsigned long int                                      element_index,
               PFuncBase<std::vector<PRISMS::Coordinate<2>>, double> *bfunc,
               const PRISMS::Coordinate<2>                           &node_coord,
               const PRISMS::Coordinate<2>                           &dim,
               int                                                    element_node_index)
      : Interpolator<Coordinate, 2>(node_index, element_index, bfunc)
    {
      _var.resize(4);

      _var[1][0] = node_coord[0];
      _var[1][1] = node_coord[1];

      _var[2][0] = dim[0];
      _var[2][1] = dim[1];

      if (element_node_index == 0)
        {
          _var[3][0] = 1.0;
          _var[3][1] = 1.0;
        }
      else if (element_node_index == 1)
        {
          _var[3][0] = -1.0;
          _var[3][1] = 1.0;
        }
      else if (element_node_index == 2)
        {
          _var[3][0] = -1.0;
          _var[3][1] = -1.0;
        }
      else if (element_node_index == 3)
        {
          _var[3][0] = 1.0;
          _var[3][1] = -1.0;
        }
    }

    [[nodiscard]] PRISMS::Coordinate<2>
    min() const override
    {
      PRISMS::Coordinate<2> coord = _var[1];

      if (_var[3][0] == -1.0)
        {
          coord[0] -= _var[2][0];
        }
      if (_var[3][1] == -1.0)
        {
          coord[1] -= _var[2][1];
        }

      return coord;
    }

    [[nodiscard]] PRISMS::Coordinate<2>
    max() const override
    {
      PRISMS::Coordinate<2> coord = _var[1];

      if (_var[3][0] == 1.0)
        {
          coord[0] += _var[2][0];
        }
      if (_var[3][1] == 1.0)
        {
          coord[1] += _var[2][1];
        }

      return coord;
    }

    bool
    is_in_range(const Coordinate &coord) override
    {
      _var[0][0] = coord[0];
      _var[0][1] = coord[1];

      for (int i = 0; i < 2; i++)
        {
          double e = _var[3][i] * (_var[0][i] - _var[1][i]) / _var[2][i];
          if (e < 0.0 || e >= 1.0)
            {
              return false;
            }
        }

      return true;
    }

    // for the following, you are expected to KNOW that the coord is_in_range!!!

    double
    operator()(const Coordinate &coord) override
    {
      _var[0][0] = coord[0];
      _var[0][1] = coord[1];
      return (*this->_bfunc)(_var);
    }
  };

} // namespace PRISMS

#endif