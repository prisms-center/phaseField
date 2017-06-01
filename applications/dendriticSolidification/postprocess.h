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
#define pp_variable_name {"theta"}
#define pp_variable_type {"SCALAR"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define pp_need_val {true,true,true}
#define pp_need_grad {true,false,true}
#define pp_need_hess  {false,false,false}

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

	// The concentration and its derivatives (names here should match those in the macros above)
	scalarvalueType c = modelVariablesList[0].scalarValue;
	dealii::Tensor<1, dim, dealii::VectorizedArray<double> > cx = modelVariablesList[0].scalarGrad;

	// The chemical potential and its derivatives (names here should match those in the macros above)
	scalarvalueType mu = modelVariablesList[1].scalarValue;

	// The order parameter and its derivatives (names here should match those in the macros above)
	scalarvalueType n = modelVariablesList[2].scalarValue;
	dealii::Tensor<1, dim, dealii::VectorizedArray<double> > nx = modelVariablesList[2].scalarGrad;

	scalarvalueType theta;

	for (unsigned i=0; i<c.n_array_elements;i++){
//		if (std::abs(nx[0][i]) < 1.0e-10) {
//			theta[i] = 1.0e-10;
//		}else{
//			theta[i] = std::atan(nx[1][i]/nx[0][i]);
//		}
		theta[i] = std::atan2(nx[1][i],nx[0][i]);
	}

	modelResidualsList[0].scalarValueResidual = WV;

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
