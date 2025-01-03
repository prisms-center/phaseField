#ifndef PField_HH
#define PField_HH

#include <field_input/IntegrationTools/pfield/Mesh.hh>
#include <field_input/IntegrationTools/pfunction/PFuncBase.hh>

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

    FieldType _val;

    // ----------------------------------------------------------
    // Constructors
    // PField();

    PField(const std::string            &name,
           Mesh<Coordinate, DIM>        &mesh,
           const std::vector<FieldType> &field,
           const FieldType              &zero)
      : PFuncBase<Coordinate, FieldType>(name)
      , _mesh(&mesh)
      , _field(field)
      , _zero(zero)

    {
      int max = mesh.max_bin_size();
      _bfunc.resize(max);
      _node_index.resize(max);
    }

    // ----------------------------------------------------------
    // Use these functions if you want to evaluate a single value
    FieldType
    operator()(const Coordinate &coord) override;

    // ----------------------------------------------------------
    // Use these functions to evaluate several values, then use 'get' methods to access
    // results
    void
    eval(const Coordinate &coord) override;

    FieldType
    operator()() const override;

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
    // get evaluated basis functions
    (*_mesh).basis_functions(coord, _bfunc, _node_index, _Nbfunc);

    // sum them
    _val = _zero;
    for (int i = 0; i < _Nbfunc; i++)
      {
        _val += _bfunc[i] * _field[_node_index[i]];
      }

    return _val;
  }

  template <class Coordinate, class FieldType, int DIM>
  void
  PField<Coordinate, FieldType, DIM>::eval(const Coordinate &coord)
  {
    (*this)(coord);
  }

  template <class Coordinate, class FieldType, int DIM>
  FieldType
  PField<Coordinate, FieldType, DIM>::operator()() const
  {
    return _val;
  }

} // namespace PRISMS

#endif
