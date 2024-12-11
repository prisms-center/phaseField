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
