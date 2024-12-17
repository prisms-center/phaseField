#ifndef Interpolator_HH
#define Interpolator_HH

#include <field_input/IntegrationTools/pfield/Coordinate.hh>
#include <field_input/IntegrationTools/pfunction/PFuncBase.hh>

namespace PRISMS
{

  /// A base class for interpolating functions
  ///
  template <class Coordinate, int DIM>
  class Interpolator
  {
  protected:
    unsigned long int _node;    // index of nodal value or control point
    unsigned long int _element; // index of element
    PFuncBase<std::vector<PRISMS::Coordinate<DIM>>, double> *_bfunc; // basis function to
                                                                     // evaluate

  public:
    using size_type =
      typename PFuncBase<std::vector<PRISMS::Coordinate<DIM>>, double>::size_type;

    Interpolator(unsigned long int                                        node,
                 unsigned long int                                        element,
                 PFuncBase<std::vector<PRISMS::Coordinate<DIM>>, double> *bfunc)
      : _node(node)
      , _element(element)
      , _bfunc(bfunc) {};

    virtual ~Interpolator() {};

    unsigned long int
    node()
    {
      return _node;
    }

    unsigned long int
    element()
    {
      return _element;
    }

    virtual PRISMS::Coordinate<DIM>
    min() const
    {
      undefined("void min(Coordinate &coord) const");
      return PRISMS::Coordinate<DIM>();
    }

    virtual PRISMS::Coordinate<DIM>
    max() const
    {
      undefined("void max(Coordinate &coord) const");
      return PRISMS::Coordinate<DIM>();
    }

    virtual bool
    is_in_range(const Coordinate &coord)
    {
      undefined("bool is_in_range(Coordinate coord) const");
      return false;
    }

    virtual double
    operator()(const Coordinate &coord)
    {
      undefined("double operator()(Coordinate coord)");
      return double();
    }

  private:
    void
    undefined(std::string fname) const
    {
      std::cout << "Error in Interpolator." << std::endl;
      std::cout << "   The member function '" << fname << "' has not been defined."
                << std::endl;
      exit(1);
    }
  };

} // namespace PRISMS

#endif
