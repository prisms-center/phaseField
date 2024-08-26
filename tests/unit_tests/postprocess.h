/*
 * postprocess.h
 *
 *  Created on: Mar 29, 2017
 *      Author: stephendewitt
 */

#ifndef APPLICATIONS_UNITTESTS_POSTPROCESS_H_
#define APPLICATIONS_UNITTESTS_POSTPROCESS_H_

// =================================================================================

template <int dim>
void
postProcessedFields(const std::vector<modelVariable<dim>> &modelVariablesList,
                    std::vector<modelResidual<dim>>       &modelResidualsList,
                    const dealii::Point<dim, dealii::VectorizedArray<double>> q_point_loc,
                    const list_of_CIJ<dim> material_moduli)
{}

// =================================================================================

// Template instantiations
template void
postProcessedFields<2>(
  const std::vector<modelVariable<2>>                    &modelVariablesList,
  std::vector<modelResidual<2>>                          &modelResidualsList,
  const dealii::Point<2, dealii::VectorizedArray<double>> q_point_loc,
  const list_of_CIJ<2>                                    material_moduli);

template void
postProcessedFields<3>(
  const std::vector<modelVariable<3>>                    &modelVariablesList,
  std::vector<modelResidual<3>>                          &modelResidualsList,
  const dealii::Point<3, dealii::VectorizedArray<double>> q_point_loc,
  const list_of_CIJ<3>                                    material_moduli);

#endif /* APPLICATIONS_UNITTESTS_POSTPROCESS_H_ */
