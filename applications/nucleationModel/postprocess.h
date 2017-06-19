/*
 * postprocess.h
 *
 *  Created on: Mar 29, 2017
 *      Author: stephendewitt
 */

#ifndef APPLICATIONS_ALLENCAHN_POSTPROCESS_H_
#define APPLICATIONS_ALLENCAHN_POSTPROCESS_H_

#include "../../include/list_of_CIJ.h"

// =================================================================================
// Define the variables for postprocessing
// =================================================================================

//// The names of the variables, whether they are scalars or vectors and whether the
//// governing eqn for the variable is parabolic or elliptic
//#define pp_variable_name {}
//#define pp_variable_type {}
//
//// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
//#define pp_need_val {}
//#define pp_need_grad {}
//#define pp_need_hess  {}
//
//// Flags for whether the residual equation has a term multiplied by the test function
//// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
//#define pp_need_val_residual {}
//#define pp_need_grad_residual {}

// =================================================================================

template <int dim>
void postProcessedFields(const std::vector<modelVariable<dim> > & modelVariablesList,
												std::vector<modelResidual<dim> > & modelResidualsList,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc,
												const list_of_CIJ<dim> material_moduli) {

// The order parameter and its derivatives (names here should match those in the macros above)




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
