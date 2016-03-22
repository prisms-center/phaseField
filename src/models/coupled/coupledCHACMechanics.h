//Matrix Free implementation of coupled Cahn-Hilliard, Allen-Cahn and Mechanics formulation 
#ifndef CHACMECHANICS_H
#define CHACMECHANICS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

//material models
#include "../mechanics/computeStress.h"

template <int dim>
class CoupledCHACMechanicsProblem: public MatrixFreePDE<dim>
{
 public: 
  CoupledCHACMechanicsProblem();

 private:
  //elasticity matrix
  Table<2, double> CIJ;

  //RHS implementation for explicit solve
  void getRHS(const MatrixFree<dim,double> &data, 
	      std::vector<vectorType*> &dst, 
	      const std::vector<vectorType*> &src,
	      const std::pair<unsigned int,unsigned int> &cell_range) const;
    
  //LHS implementation for implicit solve 
  void  getLHS(const MatrixFree<dim,double> &data, 
	       vectorType &dst, 
	       const vectorType &src,
	       const std::pair<unsigned int,unsigned int> &cell_range) const;

  //method to apply initial conditions
  void applyInitialConditions();
 
  //methods to apply dirichlet BC's on displacement
  void applyDirichletBCs();
};

//constructor
template <int dim>
CoupledCHACMechanicsProblem<dim>::CoupledCHACMechanicsProblem(): MatrixFreePDE<dim>(),
  CIJ(2*dim-1+dim/3,2*dim-1+dim/3)
{
  //initialize elasticity matrix
#if defined(MaterialModelV) && defined(MaterialConstantsV)
  double materialConstants[]=MaterialConstantsV;
  getCIJMatrix<dim>(MaterialModelV, materialConstants, CIJ, this->pcout);
#else
#error Compile ERROR: missing material property variable: MaterialModelV, MaterialConstantsV
#endif
}

template <int dim>
void  CoupledCHACMechanicsProblem<dim>::getRHS(const MatrixFree<dim,double> &data, 
					       std::vector<vectorType*> &dst, 
					       const std::vector<vectorType*> &src,
					       const std::pair<unsigned int,unsigned int> &cell_range) const{
  
  //initialize fields
  typeScalar cVals(data, 0), n1Vals(data,1), n2Vals(data,2), n3Vals(data,3);
  typeVector uVals(data, 4);

  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //initialize c field
    cVals.reinit(cell); cVals.read_dof_values_plain(*src[0]); cVals.evaluate(true, true, false);

    //initialize n fields
    n1Vals.reinit(cell); n1Vals.read_dof_values_plain(*src[1]); n1Vals.evaluate(true, true, false);
    n2Vals.reinit(cell); n2Vals.read_dof_values_plain(*src[2]); n2Vals.evaluate(true, true, false);
    n3Vals.reinit(cell); n3Vals.read_dof_values_plain(*src[3]); n3Vals.evaluate(true, true, false);

    //initialize u field 
    uVals.reinit(cell); uVals.read_dof_values_plain(*src[4]);

    if (c_dependent_misfit == true){
    	uVals.evaluate(false, true, true);
    }
    else{
    	uVals.evaluate(false, true, false);
    }


    //loop over quadrature points
    for (unsigned int q=0; q<cVals.n_q_points; ++q){
      //c
      scalarvalueType c = cVals.get_value(q);
      scalargradType cx = cVals.get_gradient(q);

      //n1
      scalarvalueType n1 = n1Vals.get_value(q);
      scalargradType n1x = n1Vals.get_gradient(q);

      //n2
      scalarvalueType n2 = n2Vals.get_value(q);
      scalargradType n2x = n2Vals.get_gradient(q);

      //n3 
      scalarvalueType n3 = n3Vals.get_value(q);
      scalargradType n3x = n3Vals.get_gradient(q);
      
      //u
      vectorgradType ux = uVals.get_gradient(q);
      vectorgradType Rux;
      vectorhessType uxx;

      if (c_dependent_misfit == true){
    	  uxx = uVals.get_hessian(q);
      }
      
      // Calculate the stress-free transformation strain and its derivatives at the quadrature point
      dealii::VectorizedArray<double> sfts1[dim][dim], sfts1c[dim][dim], sfts1cc[dim][dim], sfts2[dim][dim], sfts2c[dim][dim], sfts2cc[dim][dim], sfts3[dim][dim], sfts3c[dim][dim], sfts3cc[dim][dim];

      for (unsigned int i=0; i<dim; i++){
    	  for (unsigned int j=0; j<dim; j++){
    		  if (c_dependent_misfit == true){
    			  // Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c + b_p
    			  sfts1[i][j] = constV(a1[i][j])*c + constV(b1[i][j]);
    			  sfts1c[i][j] = constV(a1[i][j]);
    			  sfts1cc[i][j] = constV(0.0);

    			  // Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c + b_p
    			  sfts2[i][j] = constV(a2[i][j])*c + constV(b2[i][j]);
    			  sfts2c[i][j] = constV(a2[i][j]);
    			  sfts2cc[i][j] = constV(0.0);

    			  // Polynomial fits for the stress-free transformation strains, of the form: sfts = a_p * c + b_p
    			  sfts3[i][j] = constV(a3[i][j])*c + constV(b3[i][j]);
    			  sfts3c[i][j] = constV(a3[i][j]);
    			  sfts3cc[i][j] = constV(0.0);
    		  }
    		  else{
    			  sfts1[i][j] = constV(sf1Strain[i][j]);
    			  sfts1c[i][j] = constV(0.0);
    			  sfts1cc[i][j] = constV(0.0);

    			  sfts2[i][j] = constV(sf2Strain[i][j]);
    			  sfts2c[i][j] = constV(0.0);
    			  sfts2cc[i][j] = constV(0.0);

    			  sfts3[i][j] = constV(sf3Strain[i][j]);
    			  sfts3c[i][j] = constV(0.0);
    			  sfts3cc[i][j] = constV(0.0);
    		  }
    	  }
      }


      //compute E2=(E-E0)
      dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];
      for (unsigned int i=0; i<dim; i++){
    	  for (unsigned int j=0; j<dim; j++){
    		  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-(sfts1[i][j]*h1V+sfts2[i][j]*h2V+sfts3[i][j]*h3V);
    	  }
      }
      
      //compute stress
      //S=C*(E-E0)
      //Table<2, double> CIJ;
		if (n_dependent_stiffness == true){
			for (unsigned int i=0; i<2*dim-1+dim/3; i++){
				for (unsigned int j=0; j<2*dim-1+dim/3; j++){
					//CIJ(i,j) = CIJ_alpha(i,j); //CIJ_alpha(i,j)*(1.0-h1V[0]-h2V[0]-h3V[0]) + CIJ_beta(i,j)*(h1V[0]+h2V[0]+h3V[0]);
			  }
			}
		}
		else{
			for (unsigned int i=0; i<2*dim-1+dim/3; i++){
				for (unsigned int j=0; j<2*dim-1+dim/3; j++){
					//CIJ(i,j) = CIJ_homo(i,j);
				}
			}
		}
      computeStress<dim>(CIJ, E2, S);
      
      //fill residual corresponding to mechanics
      //R=-C*(E-E0)
      for (unsigned int i=0; i<dim; i++){
    	  for (unsigned int j=0; j<dim; j++){
    		  Rux[i][j] -= S[i][j];
    	  }
      }
      
      // Compute one of the stress terms in the order parameter chemical potential, nDependentMisfitACp = C*(E-E0)*(E0_p*Hn)
      dealii::VectorizedArray<double> nDependentMisfitAC1=make_vectorized_array(0.0);
      dealii::VectorizedArray<double> nDependentMisfitAC2=make_vectorized_array(0.0);
      dealii::VectorizedArray<double> nDependentMisfitAC3=make_vectorized_array(0.0);

      for (unsigned int i=0; i<dim; i++){
    	  for (unsigned int j=0; j<dim; j++){
    		  nDependentMisfitAC1+=S[i][j]*(sfts1[i][j]);
    		  nDependentMisfitAC2+=S[i][j]*(sfts2[i][j]);
    		  nDependentMisfitAC3+=S[i][j]*(sfts3[i][j]);
    	  }
      }

      nDependentMisfitAC1*=-hn1V;
      nDependentMisfitAC2*=-hn2V;
      nDependentMisfitAC3*=-hn3V;

      // Compute the other stress term in the order parameter chemical potential, heterMechACp = 0.5*Hn*(C_alpha-C_beta)*(E-E0)*(E-E0)
		dealii::VectorizedArray<double> heterMechAC1=constV(0.0);
		dealii::VectorizedArray<double> heterMechAC2=constV(0.0);
		dealii::VectorizedArray<double> heterMechAC3=constV(0.0);
		dealii::VectorizedArray<double> S2[dim][dim];

  //      if (n_dependent_stiffness == true){
  //    	  Table<2, double> CIJ_diff;
  //    	  for (unsigned int i=0; i<2*dim-1+dim/3; i++){
  //    		  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
  //    			  CIJ_diff(i,j) = 0.0; //CIJ_beta(i,j)- CIJ_alpha(i,j);
  //    		  }
  //    	  }
  //    	  computeStress<dim>(CIJ_diff, E2, S2);
  //    	  for (unsigned int i=0; i<dim; i++){
  //    		  for (unsigned int j=0; j<dim; j++){
  //    			  heterMechAC1 += S[i][j]*E2[i][j];
  //    		  }
  //    	  }
  //    	  // Aside from HnpV, heterMechAC1, heterMechAC2, and heterMechAC3 are equal
  //    	  heterMechAC2 = 0.5*hn2V*heterMechAC1;
  //    	  heterMechAC3 = 0.5*hn3V*heterMechAC1;
  //
  //    	  heterMechAC1 = 0.5*hn1V*heterMechAC1;
  //      }


		// compute the stress term in the gradient of the concentration chemical potential, grad_mu_el = [C*(E-E0)*E0c]x, must be a vector with length dim
		dealii::VectorizedArray<double> grad_mu_el[dim];

		for (unsigned int k=0; k<dim; k++){
		  grad_mu_el[k] = constV(0.0);
		}

		if (c_dependent_misfit == true){
		  dealii::VectorizedArray<double> E3[dim][dim], S3[dim][dim];

		  for (unsigned int i=0; i<dim; i++){
			  for (unsigned int j=0; j<dim; j++){
				  E3[i][j] =  -( sfts1c[i][j]*h1V + sfts2c[i][j]*h2V + sfts3c[i][j]*h3V);
			  }
		  }

		  computeStress<dim>(CIJ, E3, S3);

		  for (unsigned int i=0; i<dim; i++){
			  for (unsigned int j=0; j<dim; j++){
				  for (unsigned int k=0; k<dim; k++){
					  grad_mu_el[k]+=S3[i][j] * (constV(0.5)*(uxx[i][j][k]+uxx[j][i][k]) - (sfts1c[i][j]*h1V + sfts2c[i][j]*h2V + sfts3c[i][j]*h3V)*cx[k]
										  - (sfts1[i][j]*hn1V*n1x[k] + sfts2[i][j]*hn2V*n2x[k] + sfts3[i][j]*hn3V*n3x[k]));
				  }
			  }
		  }

		  for (unsigned int i=0; i<dim; i++){
			  for (unsigned int j=0; j<dim; j++){
				  for (unsigned int k=0; k<dim; k++){
					  grad_mu_el[k]+= - S[i][j] * (sfts1c[i][j]*hn1V*n1x[k] + sfts2c[i][j]*hn2V*n2x[k] + sfts3c[i][j]*hn3V*n3x[k]
										  + (sfts1cc[i][j]*h1V + sfts2cc[i][j]*h2V + sfts3cc[i][j]*h3V)*cx[k]);
				  }
			  }
		  }

  //    	  for (unsigned int i=0; i<dim; i++){
  //    		  for (unsigned int j=0; j<dim; j++){
  //    			  for (unsigned int k=0; k<dim; k++){
  //    				  grad_mu_el[k]+= - S2[i][j] * (sfts1c[i][j]*hn1V*n1x[k] + sfts2c[i][j]*hn2V*n2x[k] + sfts3c[i][j]*hn3V*n3x[k]);
  //    			  }
  //    		  }
  //    	  }
		}


      //compute K*nx
      scalargradType Knx1, Knx2, Knx3;
      for (unsigned int a=0; a<dim; a++) {
	Knx1[a]=0.0;
	Knx2[a]=0.0;
	Knx3[a]=0.0;
	for (unsigned int b=0; b<dim; b++){
	  Knx1[a]+=Kn1[a][b]*n1x[b];
	  Knx2[a]+=Kn2[a][b]*n2x[b];
	  Knx3[a]+=Kn3[a][b]*n3x[b];
	}
      }
  
      //submit values
      cVals.submit_value(rcV,q); cVals.submit_gradient(rcxV,q);
      n1Vals.submit_value(rn1V,q); n1Vals.submit_gradient(rn1xV,q);
      n2Vals.submit_value(rn2V,q); n2Vals.submit_gradient(rn2xV,q);
      n3Vals.submit_value(rn3V,q); n3Vals.submit_gradient(rn3xV,q);
      uVals.submit_gradient(Rux,q);
    }
    
    //integrate values
    cVals.integrate(true, true);  cVals.distribute_local_to_global(*dst[0]);
    n1Vals.integrate(true, true); n1Vals.distribute_local_to_global(*dst[1]);
    n2Vals.integrate(true, true); n2Vals.distribute_local_to_global(*dst[2]);
    n3Vals.integrate(true, true); n3Vals.distribute_local_to_global(*dst[3]);
    uVals.integrate(false, true); uVals.distribute_local_to_global(*dst[4]);
  }
}

template <int dim>
void  CoupledCHACMechanicsProblem<dim>::getLHS(const MatrixFree<dim,double> &data, 
					       vectorType &dst, 
					       const vectorType &src,
					       const std::pair<unsigned int,unsigned int> &cell_range) const{
  typeVector uVals(data, 4);
  
  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //initialize u field 
    uVals.reinit(cell); uVals.read_dof_values_plain(src); uVals.evaluate(false, true, false);

    //initialize n fields
    //n1Vals.reinit(cell); n1Vals.read_dof_values_plain(*src[1]); n1Vals.evaluate(true, false, false);
    //n2Vals.reinit(cell); n2Vals.read_dof_values_plain(*src[2]); n2Vals.evaluate(true, false, false);
    //n3Vals.reinit(cell); n3Vals.read_dof_values_plain(*src[3]); n3Vals.evaluate(true, false, false);

    //loop over quadrature points
    for (unsigned int q=0; q<uVals.n_q_points; ++q){
      //u
      vectorgradType ux = uVals.get_gradient(q);
      vectorgradType Rux;

      // n
      //scalarvalueType n1 = n1Vals.get_value(q);
      //scalarvalueType n2 = n2Vals.get_value(q);
      //scalarvalueType n3 = n3Vals.get_value(q);

      //compute strain tensor
      dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
      for (unsigned int i=0; i<dim; i++){
    	  for (unsigned int j=0; j<dim; j++){
    		  E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]);
    	  }
      }
    
      // Compute stress tensor
      //Table<2, double> CIJ;
      if (n_dependent_stiffness == true){
          	  for (unsigned int i=0; i<2*dim-1+dim/3; i++){
          		  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
          			  //CIJ(i,j) = CIJ_alpha(i,j)*(1.0-h1V-h2V-h3V) + CIJ_beta(i,j)*(h1V+h2V+h3V);
          		  }
          	  }
      }
      else{
    	  for (unsigned int i=0; i<2*dim-1+dim/3; i++){
    		  for (unsigned int j=0; j<2*dim-1+dim/3; j++){
    			  //CIJ(i,j) = CIJ_homo(i,j);
    		  }
    	  }
      }
      computeStress<dim>(CIJ, E, S);
      
      //compute residual
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  Rux[i][j] = S[i][j]; 
	}
      }
      
      //submit residual value
      uVals.submit_gradient(Rux,q);
    }
    
    //integrate
    uVals.integrate(false, true); uVals.distribute_local_to_global(dst);
  }
}

#endif
