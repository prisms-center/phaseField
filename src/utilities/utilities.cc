// utility functions for the MatrixFreePDE class

#include <core/matrixFreePDE.h>

// return index of given field name if exists, else throw error
template <int dim, int degree>
unsigned int
MatrixFreePDE<dim, degree>::getFieldIndex(std::string _name)
{
  for (const auto &field : fields)
    {
      if (field.name.compare(_name) == 0)
        {
          return field.index;
        }
    }
  pcout << "\nutilities.h: field '" << _name.c_str() << "' not initialized\n";
  exit(-1);
}

template class MatrixFreePDE<2, 1>;
template class MatrixFreePDE<3, 1>;

template class MatrixFreePDE<2, 2>;
template class MatrixFreePDE<3, 2>;

template class MatrixFreePDE<3, 3>;
template class MatrixFreePDE<2, 3>;

template class MatrixFreePDE<3, 4>;
template class MatrixFreePDE<2, 4>;

template class MatrixFreePDE<3, 5>;
template class MatrixFreePDE<2, 5>;

template class MatrixFreePDE<3, 6>;
template class MatrixFreePDE<2, 6>;