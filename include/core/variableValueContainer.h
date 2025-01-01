#ifndef VARIBLEVALUECONTAINER_H
#define VARIBLEVALUECONTAINER_H

#include <deal.II/base/exceptions.h>

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
  operator()(unsigned int global_variable_index) const
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
    AssertThrow(false,
                dealii::ExcMessage(
                  "PRISMS-PF Error: Attempted access of a variable value that was not "
                  "marked as needed in 'parameters.prm'. Double-check the indices in "
                  "user functions where a variable value is requested."));
  };

private:
  std::vector<unsigned int> var_index;
  std::vector<double>       value;
  unsigned int              num_entries;
};

#endif
