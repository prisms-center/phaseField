#ifndef FIELDS_H
#define FIELDS_H

#include <deal.II/base/conditional_ostream.h>

#include <core/varTypeEnums.h>
#include <stdexcept>
#include <string>
#include <utility>

/**
 * \brief Field class that handles the attributes of each field.
 */
template <int dim>
class Field
{
public:
  /**
   * \brief Constructor.
   */
  Field(fieldType _type, PDEType _pdetype, std::string _name);

  fieldType    type;
  PDEType      pdetype;
  std::string  name;
  unsigned int index;
  unsigned int startIndex;
  unsigned int numComponents;
  bool         hasDirichletBCs;
  bool         hasnonuniformDirichletBCs;
  bool         hasNeumannBCs;

private:
  /**
   * \brief Total global fields.
   */
  static unsigned int fieldCount;

  /**
   * \brief Total global indices
   */
  static unsigned int indexCount;

  /**
   * \brief Validate the PDEtype and fieldType passed into the constructor.
   */
  void
  validate_enum_types();

  /**
   * \brief Initialize the number of components based on whether the field is a scalar or
   * vector.
   */
  void
  init_components();

  /**
   * \brief Initialize the boundary condition conditionals to false.
   */
  void
  init_BCs();
};

// initialize static variables
template <int dim>
unsigned int Field<dim>::fieldCount = 0;

template <int dim>
unsigned int Field<dim>::indexCount = 0;

// constructor
template <int dim>
Field<dim>::Field(fieldType _type, PDEType _pdetype, std::string _name)
  : type(_type)
  , pdetype(_pdetype)
  , name(std::move(_name))
{
  validate_enum_types();

  // Set the index to the current field count and increment field count as new field is
  // being created
  index      = fieldCount++;
  startIndex = indexCount;

  init_components();

  init_BCs();
}

template <int dim>
void
Field<dim>::validate_enum_types()
{
  // Validate fieldType
  switch (type)
    {
      case fieldType::SCALAR:
      case fieldType::VECTOR:
        break;
      default:
        throw std::invalid_argument("Unknown field type in Field constructor.");
    }

  switch (pdetype)
    {
      case PDEType::EXPLICIT_TIME_DEPENDENT:
      case PDEType::AUXILIARY:
      case PDEType::IMPLICIT_TIME_DEPENDENT:
      case PDEType::TIME_INDEPENDENT:
        break;
      default:
        throw std::invalid_argument("Unknown PDE type in Field constructor.");
    }
}

template <int dim>
void
Field<dim>::init_components()
{
  switch (type)
    {
      case SCALAR:
        indexCount += 1;
        numComponents = 1;
        break;
      case VECTOR:
        indexCount += dim;
        numComponents = dim;
        break;
      default:
        break;
    }
}

template <int dim>
void
Field<dim>::init_BCs()
{
  // Default assignment of BCs
  hasDirichletBCs           = false;
  hasnonuniformDirichletBCs = false;
  hasNeumannBCs             = false;
}

#endif
