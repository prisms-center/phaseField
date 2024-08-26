/*
 * buildFields.cc
 *
 *  Created on: Feb 22, 2017
 *      Author: stephendewitt
 */

// =====================================================================
// FUNCTION TO BUILD THE VECTOR OF FIELDS
// =====================================================================

#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::buildFields()
{
  // Build each of the fields in the system
  for (unsigned int i = 0; i < userInputs.number_of_variables; i++)
    {
      fields.push_back(Field<dim>(userInputs.var_type[i],
                                  userInputs.var_eq_type[i],
                                  userInputs.var_name[i]));
    }
}

#include "../../include/matrixFreePDE_template_instantiations.h"
