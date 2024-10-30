
#ifndef Quad_HH
#define Quad_HH

#include "../../pfunction/PFuncBase.hh"
#include "../../pfunction/PSimpleBase.hh"
#include "../../pfunction/PSimpleFunction.hh"
#include "../Coordinate.hh"
#include "./Interpolator.hh"

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

    [[nodiscard]] Quad_f *
    clone() const override
    {
      return new Quad_f(*this);
    }
  };

  class Quad_grad_0 : public PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double>
  {
    [[nodiscard]] double
    eval(const std::vector<PRISMS::Coordinate<2>> &var) const override
    {
      return -var[3][0] * (1.0 - var[3][1] * (var[0][1] - var[1][1]) / var[2][1]) /
             var[2][0];
    }

  public:
    Quad_grad_0()
    {
      this->_name = "Quad_grad_0";
    }

    [[nodiscard]] Quad_grad_0 *
    clone() const override
    {
      return new Quad_grad_0(*this);
    }
  };

  class Quad_grad_1 : public PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double>
  {
    [[nodiscard]] double
    eval(const std::vector<PRISMS::Coordinate<2>> &var) const override
    {
      return -var[3][1] * (1.0 - var[3][0] * (var[0][0] - var[1][0]) / var[2][0]) /
             var[2][1];
    }

  public:
    Quad_grad_1()
    {
      this->_name = "Quad_grad_1";
    }

    [[nodiscard]] Quad_grad_1 *
    clone() const override
    {
      return new Quad_grad_1(*this);
    }
  };

  class Quad_hess_0_0 : public PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double>
  {
    [[nodiscard]] double
    eval(const std::vector<PRISMS::Coordinate<2>> &var) const override
    {
      return 0.0;
    }

  public:
    Quad_hess_0_0()
    {
      this->_name = "Quad_hess_0_0";
    }

    [[nodiscard]] Quad_hess_0_0 *
    clone() const override
    {
      return new Quad_hess_0_0(*this);
    }
  };

  class Quad_hess_0_1 : public PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double>
  {
    [[nodiscard]] double
    eval(const std::vector<PRISMS::Coordinate<2>> &var) const override
    {
      return var[3][0] * var[3][1] / var[2][0] / var[2][1];
    }

  public:
    Quad_hess_0_1()
    {
      this->_name = "Quad_hess_0_1";
    }

    [[nodiscard]] Quad_hess_0_1 *
    clone() const override
    {
      return new Quad_hess_0_1(*this);
    }
  };

  class Quad_hess_1_0 : public PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double>
  {
    [[nodiscard]] double
    eval(const std::vector<PRISMS::Coordinate<2>> &var) const override
    {
      return var[3][0] * var[3][1] / var[2][0] / var[2][1];
    }

  public:
    Quad_hess_1_0()
    {
      this->_name = "Quad_hess_1_0";
    }

    [[nodiscard]] Quad_hess_1_0 *
    clone() const override
    {
      return new Quad_hess_1_0(*this);
    }
  };

  class Quad_hess_1_1 : public PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double>
  {
    [[nodiscard]] double
    eval(const std::vector<PRISMS::Coordinate<2>> &var) const override
    {
      return 0.0;
    }

  public:
    Quad_hess_1_1()
    {
      this->_name = "Quad_hess_1_1";
    }

    [[nodiscard]] Quad_hess_1_1 *
    clone() const override
    {
      return new Quad_hess_1_1(*this);
    }
  };

  class Quad : public PFuncBase<std::vector<PRISMS::Coordinate<2>>, double>
  {
    PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double>   *_val;
    PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double>  **_grad_val;
    PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double> ***_hess_val;

  public:
    typedef PFuncBase<std::vector<PRISMS::Coordinate<2>>, double>::size_type size_type;

    Quad()
      : PFuncBase<std::vector<PRISMS::Coordinate<2>>, double>()
    {
      construct();
    }

    Quad(const Quad &RHS)
      : PFuncBase<std::vector<PRISMS::Coordinate<2>>, double>(RHS)
    {
      construct();
    }

    ~Quad()
    {
      delete _val;

      delete _grad_val[0];
      delete _grad_val[1];
      delete[] _grad_val;

      delete _hess_val[0][0];
      delete _hess_val[0][1];
      delete _hess_val[1][0];
      delete _hess_val[1][1];
      delete[] _hess_val[0];
      delete[] _hess_val[1];
      delete[] _hess_val;
    }

    [[nodiscard]] Quad *
    clone() const override
    {
      return new Quad(*this);
    }

    PSimpleFunction<std::vector<PRISMS::Coordinate<2>>, double>
    simplefunction() const override
    {
      return PSimpleFunction<std::vector<PRISMS::Coordinate<2>>, double>(*_val);
    }

    PSimpleFunction<std::vector<PRISMS::Coordinate<2>>, double>
    grad_simplefunction(size_type di) const override
    {
      return PSimpleFunction<std::vector<PRISMS::Coordinate<2>>, double>(*_grad_val[di]);
    }

    PSimpleFunction<std::vector<PRISMS::Coordinate<2>>, double>
    hess_simplefunction(size_type di, size_type dj) const override
    {
      return PSimpleFunction<std::vector<PRISMS::Coordinate<2>>, double>(
        *_hess_val[di][dj]);
    }

    double
    operator()(const std::vector<PRISMS::Coordinate<2>> &var) override
    {
      return (*_val)(var);
    }

    double
    grad(const std::vector<PRISMS::Coordinate<2>> &var, size_type di) override
    {
      return (*_grad_val[di])(var);
    }

    double
    hess(const std::vector<PRISMS::Coordinate<2>> &var,
         size_type                                 di,
         size_type                                 dj) override
    {
      return (*_hess_val[di][dj])(var);
    }

    void
    eval(const std::vector<PRISMS::Coordinate<2>> &var) override
    {
      (*_val)(var);
    }

    void
    eval_grad(const std::vector<PRISMS::Coordinate<2>> &var) override
    {
      (*_grad_val[0])(var);
      (*_grad_val[1])(var);
    }

    void
    eval_hess(const std::vector<PRISMS::Coordinate<2>> &var) override
    {
      (*_hess_val[0][0])(var);
      (*_hess_val[0][1])(var);
      (*_hess_val[1][0])(var);
      (*_hess_val[1][1])(var);
    }

    double
    operator()() const override
    {
      return (*_val)();
    }

    double
    grad(size_type di) const override
    {
      return (*_grad_val[di])();
    }

    double
    hess(size_type di, size_type dj) const override
    {
      return (*_hess_val[di][dj])();
    }

  private:
    void
    construct()
    {
      this->_name = "Quad";
      this->_var_name.clear();
      this->_var_name.push_back("r");
      this->_var_name.push_back("n");
      this->_var_name.push_back("h");
      this->_var_name.push_back("s");
      this->_var_description.clear();
      this->_var_description.push_back("Coordinate to be evaluated (Cartesian)");
      this->_var_description.push_back("Coordinate of node");
      this->_var_description.push_back("Coordinate containing element dimensions");
      this->_var_description.push_back(
        "Coordinate containing +/- 1.0, depending on which corner of quad element");

      _val = new Quad_f();

      _grad_val    = new PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double> *[2];
      _grad_val[0] = new Quad_grad_0();
      _grad_val[1] = new Quad_grad_1();

      _hess_val       = new PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double> **[2];
      _hess_val[0]    = new PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double> *[2];
      _hess_val[1]    = new PSimpleBase<std::vector<PRISMS::Coordinate<2>>, double> *[2];
      _hess_val[0][0] = new Quad_hess_0_0();
      _hess_val[0][1] = new Quad_hess_0_1();
      _hess_val[1][0] = new Quad_hess_1_0();
      _hess_val[1][1] = new Quad_hess_1_1();
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
          //_var[1][0] = node_coord[0];
          //_var[1][1] = node_coord[1];

          _var[3][0] = 1.0;
          _var[3][1] = 1.0;
        }
      else if (element_node_index == 1)
        {
          //_var[1][0] = node_coord[0] + _var[2][0];
          //_var[1][1] = node_coord[1];

          _var[3][0] = -1.0;
          _var[3][1] = 1.0;
        }
      else if (element_node_index == 2)
        {
          //_var[1][0] = node_coord[0] + _var[2][0];
          //_var[1][1] = node_coord[1] + _var[2][1];

          _var[3][0] = -1.0;
          _var[3][1] = -1.0;
        }
      else if (element_node_index == 3)
        {
          //_var[1][0] = node_coord[0];
          //_var[1][1] = node_coord[1] + _var[2][1];

          _var[3][0] = 1.0;
          _var[3][1] = -1.0;
        }
    }

    [[nodiscard]] PRISMS::Coordinate<2>
    min() const override
    {
      PRISMS::Coordinate<2> coord = _var[1];

      if (_var[3][0] == -1.0)
        coord[0] -= _var[2][0];
      if (_var[3][1] == -1.0)
        coord[1] -= _var[2][1];

      return coord;
    }

    [[nodiscard]] PRISMS::Coordinate<2>
    max() const override
    {
      PRISMS::Coordinate<2> coord = _var[1];

      if (_var[3][0] == 1.0)
        coord[0] += _var[2][0];
      if (_var[3][1] == 1.0)
        coord[1] += _var[2][1];

      return coord;
    }

    bool
    is_in_range(const Coordinate &coord) override
    {
      _var[0][0] = coord[0];
      _var[0][1] = coord[1];
      double e;

      for (int i = 0; i < 2; i++)
        {
          e = _var[3][i] * (_var[0][i] - _var[1][i]) / _var[2][i];
          if (e < 0.0 || e >= 1.0)
            return false;

          // if( e == 0.0 && std::signbit(e))
          //     return false;
        }

      // std::cout << "e: " ;
      // for( int i=0; i<2; i++)
      //{
      //     e = _var[3][i]*(_var[0][i] - _var[1][i])/_var[2][i];
      //     std::cout << e << " ";
      // }
      // std::cout << std::endl;

      return true;
    }

    // for the following,
    //   you are expected to KNOW that the coord is_in_range!!!

    double
    operator()(const Coordinate &coord) override
    {
      _var[0][0] = coord[0];
      _var[0][1] = coord[1];
      return (*this->_bfunc)(_var);
    }

    double
    grad(const Coordinate &coord, size_type di) override
    {
      _var[0][0] = coord[0];
      _var[0][1] = coord[1];
      return (*this->_bfunc).grad(_var, di);
    }

    double
    hess(const Coordinate &coord, size_type di, size_type dj) override
    {
      _var[0][0] = coord[0];
      _var[0][1] = coord[1];
      return (*this->_bfunc).hess(_var, di, dj);
    }
  };

} // namespace PRISMS

#endif