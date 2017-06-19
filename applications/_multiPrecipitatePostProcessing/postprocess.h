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

// The names of the variables, whether they are scalars or vectors and whether the
// governing eqn for the variable is parabolic or elliptic
#define pp_variable_name {"f_el","mu_c","mu_n1","mu_n2","dcdt_chem","dcdt_el","dn1dt_chem","dn1dt_el"}
#define pp_variable_type {"SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR","SCALAR"}

// Flags for whether the value, gradient, and Hessian are needed in the residual eqns
#define pp_need_val {true,true,true,true,true}
#define pp_need_grad {true,true,true,true,true}
#define pp_need_hess  {false,false,false,false,true}

// Flags for whether the residual equation has a term multiplied by the test function
// (need_val_residual) and/or the gradient of the test function (need_grad_residual)
#define pp_need_val_residual {true,true,true,true,false,false,true,true}
#define pp_need_grad_residual {false,false,true,true,true,true,true,false}

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

dealii::Tensor<3, dim, dealii::VectorizedArray<double> > uxx;

bool c_dependent_misfit = false;
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		if (std::abs(sfts_linear1[i][j])>1.0e-12){
			c_dependent_misfit = true;
		}
	}
}

if (c_dependent_misfit == true){
	uxx = modelVariablesList[4].vectorHess;
}



// Calculate the derivatives of c_beta (derivatives of c_alpha aren't needed)
// Note: this section can be optimized to reduce recalculations
scalarvalueType cbn1V, cbn2V, cbn3V, cbcV, cbcn1V, cbcn2V, cbcn3V, sum_hpV, cacV;
sum_hpV = h1V+h2V+h3V;

cbcV = faccV/( (constV(1.0)-sum_hpV)*fbccV + (sum_hpV)*faccV );
cacV = fbccV/( (constV(1.0)-sum_hpV)*fbccV + (sum_hpV)*faccV );

cbn1V = hn1V * (c_alpha - c_beta) * cbcV;
cbn2V = hn2V * (c_alpha - c_beta) * cbcV;
cbn3V = hn3V * (c_alpha - c_beta) * cbcV;

cbcn1V = (faccV * (fbccV-faccV) * hn1V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic
cbcn2V = (faccV * (fbccV-faccV) * hn2V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic
cbcn3V = (faccV * (fbccV-faccV) * hn3V)/( ((1.0-sum_hpV)*fbccV + sum_hpV*faccV)*((1.0-sum_hpV)*fbccV + sum_hpV*faccV) );  // Note: this is only true if faV and fbV are quadratic


// Calculate the stress-free transformation strain and its derivatives at the quadrature point
dealii::Tensor<2, problemDIM, dealii::VectorizedArray<double> > sfts1, sfts1c, sfts1cc, sfts2, sfts2c, sfts2cc, sfts3, sfts3c, sfts3cc;

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

scalarvalueType f_el = constV(0.0);

for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  f_el += constV(0.5) * S[i][j]*E2[i][j];
  }
}

// Calculate the chemical potential for the concentration
scalarvalueType mu_c = constV(0.0);
mu_c += facV*cacV * (1.0-sum_hpV) + fbcV*cbcV * sum_hpV;
for (unsigned int i=0; i<dim; i++){
  for (unsigned int j=0; j<dim; j++){
	  mu_c -= S[i][j]*( sfts1c[i][j]*h1V + sfts2c[i][j]*h2V + sfts3c[i][j]*h3V );
  }
}

// Calculate the chemical potential for the order parameters
// Compute one of the stress terms in the order parameter chemical potential, nDependentMisfitACp = -C*(E-E0)*(E0_n)
dealii::VectorizedArray<double> nDependentMisfitAC1=constV(0.0);
dealii::VectorizedArray<double> nDependentMisfitAC2=constV(0.0);

dealii::VectorizedArray<double> E4[dim][dim]; // Intermediate variable

for (unsigned int i=0; i<dim; i++){
for (unsigned int j=0; j<dim; j++){
	E4[i][j] = constV(sfts_linear1[i][j])*h1V + constV(sfts_linear2[i][j])*h2V + constV(sfts_linear3[i][j])*h3V;

	nDependentMisfitAC1 -= S[i][j]*(cbn1V*E4[i][j] + sfts1[i][j]*hn1V);
	nDependentMisfitAC2 -= S[i][j]*(cbn2V*E4[i][j] + sfts2[i][j]*hn2V);}
}

// Compute the other stress term in the order parameter chemical potential, heterMechACp = 0.5*Hn*(C_beta-C_alpha)*(E-E0)*(E-E0)
dealii::VectorizedArray<double> heterMechAC1=constV(0.0);
dealii::VectorizedArray<double> heterMechAC2=constV(0.0);
dealii::VectorizedArray<double> S2[dim][dim];

if (n_dependent_stiffness == true){
	computeStress<dim>(material_moduli.CIJ_list[1]-material_moduli.CIJ_list[0], E2, S2);
	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			heterMechAC1 += S2[i][j]*E2[i][j];
		}
	}
	heterMechAC1 *= 0.5*hn1V;

	computeStress<dim>(material_moduli.CIJ_list[2]-material_moduli.CIJ_list[0], E2, S2);
	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			heterMechAC2 += S2[i][j]*E2[i][j];
		}
	}
	heterMechAC2 *= 0.5*hn2V;
}


scalarvalueType mu_n1, mu_n2;
mu_n1 = (fbV-faV)*hn1V - (c_beta-c_alpha)*facV*hn1V + W*fbarriern1V + nDependentMisfitAC1 + heterMechAC1;
mu_n2 = (fbV-faV)*hn2V - (c_beta-c_alpha)*facV*hn2V + W*fbarriern2V + nDependentMisfitAC2 + heterMechAC2;

//compute K*nx
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > Knx1, Knx2;
for (unsigned int a=0; a<dim; a++) {
Knx1[a]=0.0;
Knx2[a]=0.0;
for (unsigned int b=0; b<dim; b++){
	  Knx1[a]+=constV(Kn1[a][b])*n1x[b];
	  Knx2[a]+=constV(Kn2[a][b])*n2x[b];
}
}

// compute the stress term in the gradient of the concentration chemical potential, grad_mu_el = [C*(E-E0)*E0c]x, must be a vector with length dim
dealii::Tensor<1, dim, dealii::VectorizedArray<double> > grad_mu_el;

if (c_dependent_misfit == true){
	dealii::VectorizedArray<double> E3[dim][dim], S3[dim][dim]; // Intermediate variables

	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			E3[i][j] =  -( sfts1c[i][j]*h1V + sfts2c[i][j]*h2V + sfts3c[i][j]*h3V );
		}
	}

	if (n_dependent_stiffness == true){
		computeStress<dim>(CIJ_combined, E3, S3);
	}
	else{
		computeStress<dim>(material_moduli.CIJ_list[0], E3, S3);
	}

	for (unsigned int i=0; i<dim; i++){
		for (unsigned int j=0; j<dim; j++){
			for (unsigned int k=0; k<dim; k++){
				grad_mu_el[k] += S3[i][j] * (constV(0.5)*(uxx[i][j][k]+uxx[j][i][k]) + E3[i][j]*cx[k]
										- ( (sfts1[i][j]*hn1V + cbn1V*E4[i][j])*n1x[k] + (sfts2[i][j]*hn2V + cbn2V*E4[i][j])*n2x[k]
										+ (sfts3[i][j]*hn3V + cbn3V*E4[i][j])*n3x[k]) );

				grad_mu_el[k]+= - S[i][j] * ( (sfts1c[i][j]*hn1V + cbcn1V*E4[i][j])*n1x[k]
										+ (sfts2c[i][j]*hn2V+ cbcn2V*E4[i][j])*n2x[k]
										+ (sfts3c[i][j]*hn3V+ cbcn3V*E4[i][j])*n3x[k]
										+ ( sfts1cc[i][j]*h1V + sfts2cc[i][j]*h2V + sfts3cc[i][j]*h3V )*cx[k]);

				if (n_dependent_stiffness == true){
					computeStress<dim>(material_moduli.CIJ_list[1]-material_moduli.CIJ_list[0], E2, S2);
					grad_mu_el[k]+= - S2[i][j] * (cbcV*E4[i][j]*hn1V*n1x[k]);
					computeStress<dim>(material_moduli.CIJ_list[2]-material_moduli.CIJ_list[0], E2, S2);
					grad_mu_el[k]+= - S2[i][j] * (cbcV*E4[i][j]*hn2V*n2x[k]);
					computeStress<dim>(material_moduli.CIJ_list[3]-material_moduli.CIJ_list[0], E2, S2);
					grad_mu_el[k]+= - S2[i][j] * (cbcV*E4[i][j]*hn3V*n3x[k]);

				}
			}
		}
	}
}


// Residuals for the equation to evolve the order parameter (names here should match those in the macros above)
modelResidualsList[0].scalarValueResidual = f_el; //constV(0.0);
modelResidualsList[1].scalarValueResidual = mu_c; //constV(0.0);
modelResidualsList[2].scalarValueResidual = mu_n1; //constV(0.0);
modelResidualsList[3].scalarValueResidual = mu_n2; //constV(0.0);

modelResidualsList[2].scalarGradResidual = Knx1;
modelResidualsList[3].scalarGradResidual = Knx2;

modelResidualsList[4].scalarGradResidual = cx + (c_alpha-c_beta)*(hn1V*n1x + hn2V*n2x + hn3V*n3x); //constV(0.0);
modelResidualsList[5].scalarGradResidual = grad_mu_el * ((h1V+h2V+h3V)*faccV+(constV(1.0)-h1V-h2V-h3V)*fbccV)/constV(faccV*fbccV); //constV(0.0);

modelResidualsList[6].scalarValueResidual = (fbV-faV)*hn1V - (c_beta-c_alpha)*facV*hn1V + W*fbarriern1V;
modelResidualsList[6].scalarGradResidual = Knx1;

modelResidualsList[7].scalarValueResidual = nDependentMisfitAC1 + heterMechAC1;
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
