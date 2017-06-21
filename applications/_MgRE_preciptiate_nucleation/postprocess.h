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

dealii::VectorizedArray<double> sum_hpV = h1V+h2V+h3V;
dealii::VectorizedArray<double> c_alpha = ((B2*c+0.5*(B1-A1)*sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV));
dealii::VectorizedArray<double> c_beta  = ((A2*c+0.5*(A1-B1)*(1.0-sum_hpV))/(A2*(sum_hpV)+B2*(1.0-sum_hpV)));

// Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
// Note: this section can be optimized to reduce recalculations
dealii::VectorizedArray<double> cbn1V, cbn2V, cbn3V, cbcV, cbcn1V, cbcn2V, cbcn3V, cacV;

cbcV = faccV/( (constV(1.0)-sum_hpV)*fbccV + (sum_hpV)*faccV );
cacV = fbccV/( (constV(1.0)-sum_hpV)*fbccV + (sum_hpV)*faccV );

cbn1V = hn1V * (c_alpha - c_beta) * cbcV;
cbn2V = hn2V * (c_alpha - c_beta) * cbcV;
cbn3V = hn3V * (c_alpha - c_beta) * cbcV;

cbcn1V = (faccV * (fbccV-faccV) * hn1V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic
cbcn2V = (faccV * (fbccV-faccV) * hn2V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic
cbcn3V = (faccV * (fbccV-faccV) * hn3V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic


// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, dim, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc, sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	// Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c_beta + b_p
	sfts1[i][j] = constV(sfts_linear1[i][j])*c_beta + constV(sfts_const1[i][j]);
	sfts1c[i][j] = constV(sfts_linear1[i][j]) * cbcV;
	sfts1cc[i][j] = constV(0.0);

	sfts2[i][j] = constV(sfts_linear2[i][j])*c_beta + constV(sfts_const2[i][j]);
	sfts2c[i][j] = constV(sfts_linear2[i][j]) * cbcV;
	sfts2cc[i][j] = constV(0.0);

	sfts3[i][j] = constV(sfts_linear3[i][j])*c_beta + constV(sfts_const3[i][j]);
	sfts3c[i][j] = constV(sfts_linear3[i][j]) * cbcV;
	sfts3cc[i][j] = constV(0.0);

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
		  CIJ_combined[i][j] = material_moduli.CIJ_list[0][i][j]*(constV(1.0)-sum_hpV) + material_moduli.CIJ_list[1][i][j]*(h1V) + material_moduli.CIJ_list[2][i][j]*(h2V) + material_moduli.CIJ_list[3][i][j]*(h3V);
	  }
}
computeStress<dim>(CIJ_combined, E2, S);
}
else{
computeStress<dim>(material_moduli.CIJ_list[0], E2, S);
}

dealii::VectorizedArray<double> f_el = constV(0.0);

for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  f_el += constV(0.5) * S[i][j]*E2[i][j];
  }
}

dealii::VectorizedArray<double> vm_stress;
if (dim == 3){
	vm_stress = (S[0][0]-S[1][1])*(S[0][0]-S[1][1]) + (S[1][1]-S[2][2])*(S[1][1]-S[2][2]) + (S[2][2]-S[0][0])*(S[2][2]-S[0][0]);
	vm_stress += constV(6.0)*(S[0][1]*S[0][1] + S[1][2]*S[1][2] + S[2][0]*S[2][0]);
	vm_stress *= constV(0.5);
	vm_stress = std::sqrt(vm_stress);
}
else {
	vm_stress = S[0][0]*S[0][0] - S[0][0]*S[1][1] + S[1][1]*S[1][1] + constV(3.0)*S[0][1]*S[0][1];
	vm_stress = std::sqrt(vm_stress);
}


// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
modelResidualsList[0].scalarValueResidual = f_el;
modelResidualsList[1].scalarValueResidual = vm_stress;


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
