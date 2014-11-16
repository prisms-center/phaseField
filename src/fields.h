//fields class 
#ifndef FIELDS_H
#define FIELDS_H

enum fieldType {SCALAR, VECTOR}
enum PDEType   {ELLIPTIC, PARABOLIC}

template<int dim>
class Field
{
 public:
  Field(fieldType _type, PDEType _pdetype, char[] _name);
  fieldType type;
  PDEType   pdetype;
  std::string name;
  unsigned int index;
  static unsigned int fieldCount;
  static unsigned int indexCount;
}

//initialize static variables
unsigned int Field::fieldCount = 0;
unsigned int Field::indexCount = 0;

//constructor
template<int dim>
void Field::Field(fieldType _type, PDEType _pdetype, char[] _name): type(_type), pdetype(_pdetype), name(_name) 
{
  //increment field count as new field is being created
  fieldCount++; 

  //assign index to this field
  switch (type){
  case SCALAR:{
    index= indexCount;
    //increment index count by one
    indexCount+=1;
    break;
  }
  case VECTOR:{
    index= indexCount;
    //increment index count by dim
    indexCount+=dim;
    break;
  }
  default:{
    std::cout << "fields.h: unknown field type\n";
    exit(-1);
  }
  }
}

#endif
