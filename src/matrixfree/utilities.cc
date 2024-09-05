// utility functions for the MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

// return index of given field name if exists, else throw error
template <int dim, int degree>
unsigned int
MatrixFreePDE<dim, degree>::getFieldIndex(std::string _name)
{
  for (const auto &field : fields)
    {
      if (field.name.compare(_name) == 0)
        return field.index;
    }
  pcout << "\nutilities.h: field '" << _name.c_str() << "' not initialized\n";
  exit(-1);
}

#include "../../include/matrixFreePDE_template_instantiations.h"
