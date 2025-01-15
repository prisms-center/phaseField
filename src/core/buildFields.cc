/*
 * buildFields.cc
 *
 *  Created on: Feb 22, 2017
 *      Author: stephendewitt
 */

// =====================================================================
// FUNCTION TO BUILD THE VECTOR OF FIELDS
// =====================================================================

#include <core/matrixFreePDE.h>

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::buildFields()
{
  // Build each of the fields in the system
  for (const auto &[index, variable] : var_attributes)
    {
      fields.push_back(Field<dim>(variable.var_type, variable.eq_type, variable.name));
    }
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