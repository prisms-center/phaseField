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
// The number of variables
#define pp_num_var 2

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define pp_variable_name {"f_el","mu_c"}
#define pp_variable_type {"SCALAR","SCALAR"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define pp_need_val {true,true,true,true,false}
#define pp_need_grad {true,true,true,true,true}
#define pp_need_hess  {false,false,false,false,false}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define pp_need_val_residual {true,true}
#define pp_need_grad_residual {false,false}

// =================================================================================

//template <int dim, int degree>
//class customPostProcessor: public PostProcessor<dim,degree>
//{
//private:

template <int dim>
void postProcessedFields(const std::vector<modelVariable<dim> > & modelVariablesList,
												std::vector<modelResidual<dim> > & modelResidualsList,
												const dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc,
												const list_of_CIJ<dim> material_moduli) {

// The order parameter and its derivatives (names here should match those in the macros above)

//c
dealii::VectorizedArray<double> c = modelVariablesList[0].scalarValue;
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > cx = modelVariablesList[0].scalarGrad;

//n1
dealii::VectorizedArray<double> n1 = modelVariablesList[1].scalarValue;
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > n1x = modelVariablesList[1].scalarGrad;

// The second order parameter and its derivatives (names here should match those in the macros above)
dealii::VectorizedArray<double> n2 = modelVariablesList[2].scalarValue;
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > n2x = modelVariablesList[2].scalarGrad;

// The third order parameter and its derivatives (names here should match those in the macros above)
dealii::VectorizedArray<double> n3 = modelVariablesList[3].scalarValue;
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > n3x = modelVariablesList[3].scalarGrad;

//u
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > ux = modelVariablesList[4].vectorGrad;


// Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
// Note: this section can be optimized to reduce recalculations
scalarvalueType cbn1V, cbn2V, cbn3V, cbcV, cbcn1V, cbcn2V, cbcn3V, sum_hpV;
sum_hpV = h1V+h2V+h3V;

cbcV = faccV/( (constV(1.0)-sum_hpV)*fbccV + (sum_hpV)*faccV );

cbn1V = hn1V * (c_alpha - c_beta) * cbcV;
cbn2V = hn2V * (c_alpha - c_beta) * cbcV;
cbn3V = hn3V * (c_alpha - c_beta) * cbcV;

cbcn1V = (faccV * (fbccV-faccV) * hn1V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic
cbcn2V = (faccV * (fbccV-faccV) * hn2V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic
cbcn3V = (faccV * (fbccV-faccV) * hn3V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic


// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc, sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc, sfts1n, sfts1cn, sfts2n, sfts2cn, sfts3n, sfts3cn;

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c_beta + b_p
	sfts1[i][j] = constV(sfts_linear1[i][j])*c_beta + constV(sfts_const1[i][j]);
	sfts1c[i][j] = constV(sfts_linear1[i][j]) * cbcV;
	sfts1cc[i][j] = constV(0.0);
	sfts1n[i][j] = constV(sfts_linear1[i][j]) * cbn1V;
	sfts1cn[i][j] = constV(sfts_linear1[i][j]) * cbcn1V;

	sfts2[i][j] = constV(sfts_linear2[i][j])*c_beta + constV(sfts_const2[i][j]);
	sfts2c[i][j] = constV(sfts_linear2[i][j]) * cbcV;
	sfts2cc[i][j] = constV(0.0);
	sfts2n[i][j] = constV(sfts_linear2[i][j]) * cbn2V;
	sfts2cn[i][j] = constV(sfts_linear2[i][j]) * cbcn2V;

	sfts3[i][j] = constV(sfts_linear3[i][j])*c_beta + constV(sfts_const3[i][j]);
	sfts3c[i][j] = constV(sfts_linear3[i][j]) * cbcV;
	sfts3cc[i][j] = constV(0.0);
	sfts3n[i][j] = constV(sfts_linear3[i][j]) * cbn3V;
	sfts3cn[i][j] = constV(sfts_linear3[i][j]) * cbcn2V;

}
}


//compute E2=(E-E0)
dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];

for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-( sfts1[i][j]*h1V + sfts2[i][j]*h2V + sfts3[i][j]*h3V);

  }
}

//compute stress
//S=C*(E-E0)
dealii::VectorizedArray<double> CIJ_combined[2*dim-1+dim/3][2*dim-1+dim/3];

if (n_dependent_stiffness == true){
  for (unsigned int i=0; i<2*dim-1+dim/3; i++){
	  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
		  CIJ_combined[i][j] = material_moduli.CIJ_list[0][i][j]*(constV(1.0)-h1strainV-h2strainV-h3strainV) + material_moduli.CIJ_list[1][i][j]*(h1strainV) + material_moduli.CIJ_list[2][i][j]*(h2strainV) + material_moduli.CIJ_list[3][i][j]*(h3strainV);
	  }
  }
  computeStress<dim>(CIJ_combined, E2, S);
}
else{
  computeStress<dim>(material_moduli.CIJ_list[0], E2, S);
}

scalarvalueType f_el = constV(0.0);

for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  f_el += constV(0.5) * S[i][j]*E2[i][j];
  }
}

scalarvalueType mu_c = constV(0.0);
mu_c += faV*cbcV * (1.0-sum_hpV) + fbV*cbcV * sum_hpV;
for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  mu_c -= S[i][j]*( sfts1c[i][j]*h1V + sfts2c[i][j]*h2V + sfts3c[i][j]*h3V );
  }
}


// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
modelResidualsList[0].scalarValueResidual = f_el; //constV(0.0);
modelResidualsList[1].scalarValueResidual = mu_c; //constV(0.0);


}
//};

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
