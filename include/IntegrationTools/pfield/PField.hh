
#ifndef PField_HH
#define PField_HH

#include "../pfunction/PFuncBase.hh"
#include "./Mesh.hh"

//_name(name),
//_var_name(var_name),
//_var_description(var_description),

namespace PRISMS
{

  /// A class for a field
  ///
  ///   Varcontainer contains variables, for instance 'std::vector<double>' for x, y, z
  ///   FieldType is the datatype for the field,
  ///       for instance 'double' for temperature or 'std::vector<double>' for vector
  ///       displacement
  ///
  ///   A field consists of a pointer to a mesh, a list field values (at mesh nodes)
  template <class Coordinate, class FieldType, int DIM>
  class PField : public PFuncBase<Coordinate, FieldType>
  {
  public:
    using size_type = typename PFuncBase<Coordinate, FieldType>::size_type;

    // pointer to a Mesh that lives in a Body
    Mesh<Coordinate, DIM> *_mesh;

    // array of field values at mesh nodes
    std::vector<FieldType> _field;

    FieldType _zero;

    FieldType                           _val;
    std::vector<FieldType>              _grad_val;
    std::vector<std::vector<FieldType>> _hess_val;

    // ----------------------------------------------------------
    // Constructors
    // PField();

    PField(const std::string              &name,
           const std::vector<std::string> &var_name,
           const std::vector<std::string> &var_description,
           Mesh<Coordinate, DIM>          &mesh,
           const std::vector<FieldType>   &field,
           const FieldType                &zero)
      : PFuncBase<Coordinate, FieldType>(name, var_name, var_description)
      , _mesh(&mesh)
      , _field(field)
      , _zero(zero)

    {
      int max = mesh.max_bin_size();
      _bfunc.resize(max);
      _node_index.resize(max);
      _grad_val.resize(DIM);
      _hess_val.resize(DIM);
      for (int i = 0; i < DIM; i++)
        _hess_val[i].resize(DIM);
    }

    // ----------------------------------------------------------
    // Clone

    PField<Coordinate, FieldType, DIM> *
    clone() const override
    {
      return new PField<Coordinate, FieldType, DIM>(*this);
    }

    // ----------------------------------------------------------
    // Use these functions if you want to evaluate a single value
    FieldType
    operator()(const Coordinate &coord) override;
    FieldType
    grad(const Coordinate &coord, size_type di) override;
    FieldType
    hess(const Coordinate &coord, size_type di, size_type dj) override;

    // ----------------------------------------------------------
    // Use these functions to evaluate several values, then use 'get' methods to access
    // results
    void
    eval(const Coordinate &coord) override;
    void
    eval_grad(const Coordinate &coord) override;
    void
    eval_hess(const Coordinate &coord) override;

    FieldType
    operator()() const override;
    FieldType
    grad(size_type di) const override;
    FieldType
    hess(size_type di, size_type dj) const override;

    // PField unique members ------------------------------------------

  private:
    // temporary vector
    std::vector<double>            _bfunc;
    std::vector<unsigned long int> _node_index;
    int                            _Nbfunc;
  };

  template <class Coordinate, class FieldType, int DIM>
  FieldType
  PField<Coordinate, FieldType, DIM>::operator()(const Coordinate &coord)
  {
    // std::cout << "begin calc" << std::endl;

    // std::cout << "coord: ";
    // for( int i=0; i<DIM; i++)
    //     std::cout << coord[i] << " ";
    // std::cout << std::endl;

    // get evaluated basis functions
    (*_mesh).basis_functions(coord, _bfunc, _node_index, _Nbfunc);

    // sum them
    _val = _zero;
    for (int i = 0; i < _Nbfunc; i++)
      {
        // std::cout << "i: " << i << "  node: " << _node_index[i] << "  bfunc: " <<
        // _bfunc[i] << "  val: " << _field[_node_index[i]] << std::endl;
        _val += _bfunc[i] * _field[_node_index[i]];
      }
    // std::cout << "val: " << _val << std::endl << std::endl;
    // std::cout << "_Nbfunc " << _Nbfunc << std::endl;
    return _val;
  }

  template <class Coordinate, class FieldType, int DIM>
  FieldType
  PField<Coordinate, FieldType, DIM>::grad(const Coordinate &coord, size_type di)
  {
    // std::cout << "begin PField::grad()" << std::endl;
    //  get evaluated basis functions
    (*_mesh).grad_basis_functions(coord, di, _bfunc, _node_index, _Nbfunc);

    // sum them
    _grad_val[di] = _zero;
    for (int i = 0; i < _Nbfunc; i++)
      _grad_val[di] += _bfunc[i] * _field[_node_index[i]];
    return _grad_val[di];
    // std::cout << "finish PField::grad()" << std::endl;
  }

  template <class Coordinate, class FieldType, int DIM>
  FieldType
  PField<Coordinate, FieldType, DIM>::hess(const Coordinate &coord,
                                           size_type         di,
                                           size_type         dj)
  {
    // get evaluated basis functions
    (*_mesh).hess_basis_functions(coord, di, dj, _bfunc, _node_index, _Nbfunc);

    // sum them
    _hess_val[di][dj] = _zero;
    for (int i = 0; i < _Nbfunc; i++)
      _hess_val[di][dj] += _bfunc[i] * _field[_node_index[i]];
    return _hess_val[di][dj];
  }

  template <class Coordinate, class FieldType, int DIM>
  void
  PField<Coordinate, FieldType, DIM>::eval(const Coordinate &coord)
  {
    (*this)(coord);
  }

  template <class Coordinate, class FieldType, int DIM>
  void
  PField<Coordinate, FieldType, DIM>::eval_grad(const Coordinate &coord)
  {
    for (int di = 0; di < DIM; di++)
      {
        (*this).grad(coord, di);
      }
  }

  template <class Coordinate, class FieldType, int DIM>
  void
  PField<Coordinate, FieldType, DIM>::eval_hess(const Coordinate &coord)
  {
    for (int di = 0; di < DIM; di++)
      {
        for (int dj = 0; dj < DIM; dj++)
          {
            (*this).hess(coord, di, dj);
          }
      }
  }

  template <class Coordinate, class FieldType, int DIM>
  FieldType
  PField<Coordinate, FieldType, DIM>::operator()() const
  {
    return _val;
  }

  template <class Coordinate, class FieldType, int DIM>
  FieldType
  PField<Coordinate, FieldType, DIM>::grad(size_type di) const
  {
    return _grad_val[di];
  }

  template <class Coordinate, class FieldType, int DIM>
  FieldType
  PField<Coordinate, FieldType, DIM>::hess(size_type di, size_type dj) const
  {
    return _hess_val[di][dj];
  }

} // namespace PRISMS

#endif
