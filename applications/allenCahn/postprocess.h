/*
 * postprocess.h
 *
 *  Created on: Mar 29, 2017
 *      Author: stephendewitt
 */

#ifndef APPLICATIONS_ALLENCAHN_POSTPROCESS_H_
#define APPLICATIONS_ALLENCAHN_POSTPROCESS_H_

// =================================================================================
// Define the variables for postprocessing
// =================================================================================
// The number of variables
#define pp_num_var 1

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define pp_variable_name {"mag_grad_n"}
#define pp_variable_type {"SCALAR"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define pp_need_val {true}
#define pp_need_grad {true}
#define pp_need_hess  {false}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define pp_need_val_residual {true}
#define pp_need_grad_residual {false}

// =================================================================================

template <int dim>
void postProcessedFields(const std::vector<modelVariable<dim> > & modelVariablesList,
												std::vector<modelResidual<dim> > & modelResidualsList,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc,
												const list_of_CIJ<dim> material_moduli) {

// The order parameter and its derivatives (names here should match those in the macros above)
dealii::VectorizedArray<double> n = modelVariablesList[0].scalarValue;
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
