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
    uVals.reinit(cell); uVals.read_dof_values_plain(*src[4]); uVals.evaluate(false, true, false);

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
      
      //compute E2=(E-E0)
      dealii::VectorizedArray<double> E2[dim][dim], S[dim][dim];
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= constV(0.5)*(ux[i][j]+ux[j][i])-(sf1Strain[i][j]*h1V+sf2Strain[i][j]*h2V+sf3Strain[i][j]*h3V);
	}
      }
      
      //compute stress
      //S=C*(E-E0)
      computeStress<dim>(CIJ, E2, S);
      
      //fill residual corresponding to mechanics
      //R=-C*(E-E0)
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  Rux[i][j] -= S[i][j]; 
	}
      }
      
      //compute C*(E-E0)*(Esf*Hn)
      VectorizedArray<double> CEE1=make_vectorized_array(0.0);
      VectorizedArray<double> CEE2=make_vectorized_array(0.0);
      VectorizedArray<double> CEE3=make_vectorized_array(0.0);
      
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  CEE1+=S[i][j]*sf1Strain[i][j];
	  CEE2+=S[i][j]*sf2Strain[i][j];
	  CEE3+=S[i][j]*sf3Strain[i][j];	
	}
      }
      CEE1*=hn1V;
      CEE2*=hn2V;
      CEE3*=hn3V;
      
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

    //loop over quadrature points
    for (unsigned int q=0; q<uVals.n_q_points; ++q){
      //u
      vectorgradType ux = uVals.get_gradient(q);
      vectorgradType Rux;

      //compute strain tensor
      dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]);
	}
      }
    
      //compute stress tensor
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
