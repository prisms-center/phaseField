//fields class 
#ifndef FIELDS_H
#define FIELDS_H
#include <deal.II/base/conditional_ostream.h>

enum fieldType {SCALAR, VECTOR};
enum PDEType   {ELLIPTIC, PARABOLIC};

template<int dim>
class Field
{
 public:
  Field(fieldType _type, PDEType _pdetype, std::string _name);
  fieldType type;
  PDEType   pdetype;
  std::string name;
  unsigned int index;
  unsigned int startIndex;
  unsigned int numComponents;

 private:
  static unsigned int fieldCount;
  static unsigned int indexCount;
};

//initialize static variables
template <int dim> unsigned int Field<dim>::fieldCount = 0;
template <int dim> unsigned int Field<dim>::indexCount = 0;

//constructor
template<int dim>
Field<dim>::Field(fieldType _type, PDEType _pdetype, std::string _name): type(_type), pdetype(_pdetype), name(_name)
{
  //increment field count as new field is being created
  index=fieldCount;
  fieldCount++; 
  startIndex= indexCount;

  //assign index to this field
  switch (type){
  case SCALAR:{
    //increment index count by one
    indexCount+=1;
    numComponents=1;
    break;
  }
  case VECTOR:{
    //increment index count by dim
    indexCount+=dim;
    numComponents=dim;
    break;
  }
  default:{
    std::cout << "fields.h: unknown field type\n";
    exit(-1);
  }
  }
}

#endif

