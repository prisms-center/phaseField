/*
 * postprocess.h
 *
 *  Created on: Mar 29, 2017
 *      Author: stephendewitt
 */

#ifndef APPLICATIONS_ALLENCAHN_POSTPROCESS_H_
#define APPLICATIONS_ALLENCAHN_POSTPROCESS_H_

template <int dim>
void postProcessedFields(const std::vector<modelVariable<dim> > & modelVariablesList,
												std::vector<modelResidual<dim> > & modelResidualsList,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc,
												const list_of_CIJ<dim> material_moduli) {

// The order parameter and its derivatives (names here should match those in the macros above)
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > nx = modelVariablesList[0].scalarGrad;

dealii::Tensor<1, dim, dealii::VectorizedArray<double> > pp_field;
pp_field[0] = nx[0];
pp_field[1] = nx[1];


// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
modelResidualsList[0].scalarValueResidual = std::sqrt(pp_field[0]*pp_field[0]+pp_field[1]*pp_field[1]); //constV(0.0);


}

// =================================================================================

// Template instantiations
template void postProcessedFields<2>(const std::vector<modelVariable<2> > & modelVariablesList,
		std::vector<modelResidual<2> > & modelResidualsList,
		const dealii::Point<2, dealii::VectorizedArray<double> > q_point_loc,
		const list_of_CIJ<2> material_moduli);

template void postProcessedFields<3>(const std::vector<modelVariable<3> > & modelVariablesList,
		std::vector<modelResidual<3> > & modelResidualsList,
		const dealii::Point<3, dealii::VectorizedArray<double> > q_point_loc,
		const list_of_CIJ<3> material_moduli);

#endif /* APPLICATIONS_ALLENCAHN_POSTPROCESS_H_ */
