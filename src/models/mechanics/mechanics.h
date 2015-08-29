//Matrix Free implementation of infinitesimal strain mechanics
#ifndef MECHANICS_MECHANICS_H
#define MECHANICS_MECHANICS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

//material models
#include "computeStress.h"

template <int dim>
class MechanicsProblem: public MatrixFreePDE<dim>
{
 public: 
  MechanicsProblem();

 private:
  //elasticity matrix
  Table<2, double> CIJ;
  
  //RHS implementation for implicit/explicit solve
  void getRHS(std::map<std::string, typeScalar*>  valsScalar, \
	      std::map<std::string, typeVector*>  valsVector, \
	      unsigned int q) const;
  
  //LHS implementation for implicit solve 
  void getLHS(typeVector& vals, unsigned int q) const;  
  
  //methods to apply dirichlet BC's
  void markBoundaries();
  void applyDirichletBCs();
};

//constructor
template <int dim>
MechanicsProblem<dim>::MechanicsProblem(): MatrixFreePDE<dim>(), 
  CIJ(2*dim-1+dim/3,2*dim-1+dim/3)
{
  //check if all required parameters correctly specified
#if numFields!=problemDIM
#error Compile ERROR: numFields!=problemDIM. Number of fields in Mechanics problem should be equal to the problem dimension
#endif
#ifndef MaterialModelV
  this->pcout << "\nERROR: missing material property variable: MaterialModelV. exiting...\n\n"; exit(-1);
#endif
#ifndef MaterialConstantsV
  this->pcout << "\nERROR: missing material property variable: MaterialConstantsV. exiting...\n\n"; exit(-1);
#endif

  //initialize elasticity matrix
#if defined(MaterialModelV) && defined(MaterialConstantsV)
  double materialConstants[]=MaterialConstantsV;
  getCIJMatrix<dim>(MaterialModelV, materialConstants, CIJ, this->pcout);
#endif

  //
  this->getValue["u"]=false; this->getGradient["u"]=true;
  this->setValue["u"]=false; this->setGradient["u"]=true;
}

//implementation of the RHS evaluation 
template <int dim>
void  MechanicsProblem<dim>::getRHS(std::map<std::string, typeScalar*>  valsScalar, \
				    std::map<std::string, typeVector*>  valsVector, \
				    unsigned int q) const{
  vectorgradType ux = valsVector["u"]->get_gradient(q);
  vectorgradType Rux;

  //compute stress
  vectorgradType S;
  computeStress<dim>(CIJ, ux, S);
  
  //fill residual
  for (unsigned int i=0; i<dim; i++){
    for (unsigned int j=0; j<dim; j++){
      Rux[i,j] = -S[i,j];
    }
  }
  
  //compute residuals
  valsVector["u"]->submit_gradient(Rux,q);
}

template <int dim>
void  MechanicsProblem<dim>::getLHS(typeVector& vals, unsigned int q) const{
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
