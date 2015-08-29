//Matrix Free implementation of coupled Cahn-Hilliard and Mechanics formulation 
#ifndef CHMECHANICS_DIFFUSION_H
#define CHMECHANICS_DIFFUSION_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

//material models
#include "../mechanics/computeStress.h"

template <int dim>
class CoupledCahnHilliardMechanicsProblem: public MatrixFreePDE<dim>
{
 public: 
  CoupledCahnHilliardMechanicsProblem();

 private:
  //elasticity matrix
  Table<2, double> CIJ;

  //RHS implementation for explicit solve
  void getRHS(std::map<std::string, typeScalar*>  valsScalar, \
	      std::map<std::string, typeVector*>  valsVector, \
	      unsigned int q) const;

  //LHS implementation for implicit solve 
  void getLHS(typeVector& vals, unsigned int q) const;  

  //method to apply initial conditions
  void applyInitialConditions();
 
  //methods to apply dirichlet BC's on displacement
  void applyDirichletBCs();
};

//constructor
template <int dim>
CoupledCahnHilliardMechanicsProblem<dim>::CoupledCahnHilliardMechanicsProblem(): MatrixFreePDE<dim>(),
  CIJ(2*dim-1+dim/3,2*dim-1+dim/3)
{
  //check if all required parameters correctly specified
#if numFields!=(problemDIM+2)
#error Compile ERROR: numFields!=(problemDIM+2). Number of fields in coupled Cahn-Hilliard and Mechanics problem should be equal to (problemDIM+2).
#endif
#if !defined(McV) || !defined(KcV)
#error Compile ERROR: missing Cahn-Hilliard parameters. Required parameters are McV (mobility) and KcV (length scale parameter).
#endif
#if !defined(rmuV) || !defined(rmuxV) || !defined(rcV) || !defined(rcxV) 
#error Compile ERROR: missing Cahn-Hilliard residual expressions. Required expressions are rmuV, rmuxV, rcV, rcxV.
#endif
#ifndef MaterialModelV
#error Compile ERROR: missing material property variable: MaterialModelV.
#endif
#ifndef MaterialConstantsV
#error  Compile ERROR: missing material property variable: MaterialConstantsV.
#endif
#ifndef chemicalStrainV
#error  Compile ERROR: missing coupled property variable: chemicalStrainV.
#endif

  //initialize elasticity matrix
#if defined(MaterialModelV) && defined(MaterialConstantsV)
  double materialConstants[]=MaterialConstantsV;
  getCIJMatrix<dim>(MaterialModelV, materialConstants, CIJ, this->pcout);
#endif

  //"c"
  this->getValue["c"]=true; this->getGradient["c"]=true;
  this->setValue["c"]=true; this->setGradient["c"]=true;
  //"mu"
  this->getValue["mu"]=true; this->getGradient["mu"]=true;
  this->setValue["mu"]=true; this->setGradient["mu"]=true;
  //"u"
  this->getValue["u"]=false; this->getGradient["u"]=true;
  this->setValue["u"]=false; this->setGradient["u"]=true;
}


template <int dim>
void  CoupledCahnHilliardMechanicsProblem<dim>::getRHS(std::map<std::string, typeScalar*>  valsScalar, \
						       std::map<std::string, typeVector*>  valsVector, \
						       unsigned int q) const{
  //"mu"fields
  scalarvalueType mu = valsScalar["mu"]->get_value(q);
  scalargradType mux = valsScalar["mu"]->get_gradient(q);

  //"c" fields
  scalarvalueType c = valsScalar["c"]->get_value(q);
  scalargradType cx = valsScalar["c"]->get_gradient(q);

  //"u" fields
  vectorgradType ux = valsVector["u"]->get_gradient(q);
  vectorgradType Rux;

  //define identity tensor
  vectorgradType I; 
  for (unsigned int i=0; i<dim; i++) {I[i][i]+=constV(1.0);}

  //apply chemical strain
  ux+=chemicalStrainV*I;

  //compute stress
  vectorgradType S;
  computeStress<dim>(CIJ, ux, S);

  //fill residual corresponding to mechanics
  for (unsigned int i=0; i<dim; i++){
    for (unsigned int j=0; j<dim; j++){
      Rux[i,j] = -S[i,j];
    }
  }
  
  //check to ensure we are working on the intended field
  //chemical potential field
  if (this->fields[this->currentFieldIndex].name.compare("mu")==0){
    //compute residuals
    valsScalar["mu"]->submit_value(rmuV,q); 
    valsScalar["mu"]->submit_gradient(rmuxV,q);
  }
  //concentration field
  else if (this->fields[this->currentFieldIndex].name.compare("c")==0){
    //compute residuals
    valsScalar["c"]->submit_value(rcV,q);   
    valsScalar["c"]->submit_gradient(rcxV,q);
  }
  //displacement field
  else if (this->fields[this->currentFieldIndex].name.compare("u")==0){
    valsVector["u"]->submit_gradient(Rux,q);
  }
}

template <int dim>
void  CoupledCahnHilliardMechanicsProblem<dim>::getLHS(typeVector& vals, unsigned int q) const{
  //check to ensure we are working on the intended implicit field
  if (this->fields[this->currentFieldIndex].name.compare("u")==0){
    vectorgradType ux = vals.get_gradient(q);
    vectorgradType Rux;
    
    //compute stress
    vectorgradType S;
    computeStress<dim>(CIJ, ux, S);
    
    //compute residual
    for (unsigned int i=0; i<dim; i++){
      for (unsigned int j=0; j<dim; j++){
	Rux[i,j] = S[i,j];
      }
    }
    
    //submit residual value for quadrature integration and assemble
    vals.submit_gradient(Rux,q);
  }
}

#endif
