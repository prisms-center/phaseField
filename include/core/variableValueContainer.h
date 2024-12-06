// This class permits the access of a subset of indexed fields and gives an
// error if any non-allowed fields are requested
#ifndef VARIBLEVALUECONTAINER_H
#define VARIBLEVALUECONTAINER_H

class variableValueContainer
{
public:
  variableValueContainer()
    : num_entries(0) {};

  void
  set(unsigned int global_variable_index, double variable_value)
  {
    var_index.push_back(global_variable_index);
    value.push_back(variable_value);
    num_entries++;
  };

  double
  operator()(unsigned int global_variable_index)
  {
    for (unsigned int i = 0; i < num_entries; i++)
      {
        if (global_variable_index == var_index[i])
          {
            return value[i];
          }
      }

    // If this point is reached, the index isn't allowed so an error should be
    // thrown
    std::cerr << "PRISMS-PF Error: Attempted access of a variable value that was not "
                 "marked as needed in 'parameters.in'. Double-check the indices in "
                 "user functions where a variable value is requested."
              << std::endl;
    abort();
  };

private:
  std::vector<unsigned int> var_index;
  std::vector<double>       value;
  unsigned int              num_entries;
};

// void variableValueContainer::set(unsigned int global_variable_index, double
// variable_value){
//     var_index.push_back(global_variable_index);
//     value.push_back(variable_value);
//     num_entries++;
// }
//
// double variableValueContainer::operator()(unsigned int
// global_variable_index){
//     for (unsigned int i=0; i < num_entries; i++){
//         if (global_variable_index == var_index[i]){
//             return value[i];
//         }
//     }
//
//     // If this point is reached, the index isn't allowed so an error should
//     be thrown std::cerr << "PRISMS-PF Error: Attempted access of a variable
//     value that was not marked as needed in 'parameters.in'. Double-check the
//     indices in user functions where a variable value is requested." <<
//     std::endl; abort();
// }

#endif
