/*
 * postProcessedFields.h
 *
 *  Created on: Mar 29, 2017
 *      Author: stephendewitt
 */

#ifndef INCLUDE_POSTPROCESSEDFIELDS_H_
#define INCLUDE_POSTPROCESSEDFIELDS_H_

#include "list_of_CIJ.h"

template <int dim>
void
postProcessedFields(const std::vector<modelVariable<dim>> &modelVariablesList,
                    std::vector<modelResidual<dim>>       &modelResidualsList,
                    const dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc,
                    const list_of_CIJ<dim> material_moduli);

#endif /* INCLUDE_POSTPROCESSEDFIELDS_H_ */
