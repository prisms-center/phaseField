//Matrix Free implementation of infinitesimal strain mechanics
#ifndef MECHANICS_MECHANICS_H
#define MECHANICS_MECHANICS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

//material models
#include "anisotropy.cc"
#include "computeStress.cc"

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
  //initialize elasticity matrix
  double materialConstants[]=MaterialConstantsv;
  getCIJMatrix<dim>(MaterialModelv, materialConstants, CIJ, this->pcout);

  //
  this->getValue["u"]=false; this->getGradient["u"]=true;
  this->setValue["u"]=false; this->setGradient["u"]=true;
}

//implementation of the RHS evaluation 
template <int dim>
void  MechanicsProblem<dim>::getRHS(std::map<std::string, typeScalar*>  valsScalar, \
				    std::map<std::string, typeVector*>  valsVector, \
				    unsigned int q) const{
  Tensor<1, dim, gradType> ux = valsVector["u"]->get_gradient(q);
  Tensor<1, dim, gradType> Rux;
  //compute stress
  dealii::VectorizedArray<double> S[dim][dim];
  computeStress<dim>(CIJ, ux, S);
  
  //fill residual
  for (unsigned int i=0; i<dim; i++){
    for (unsigned int j=0; j<dim; j++){
      Rux[1,i][1,j] = -S[i][j];
    }
  }
  
  //compute residuals
  valsVector["u"]->submit_gradient(Rux,q);
}

template <int dim>
void  MechanicsProblem<dim>::getLHS(typeVector& vals, unsigned int q) const{
  //check to ensure we are working on the intended implicit field
  if (this->fields[this->implicitFieldIndex].name.compare("u")==0){
    Tensor<1, dim, gradType> ux = vals.get_gradient(q);
    Tensor<1, dim, gradType> Rux;
    //compute stress
    dealii::VectorizedArray<double> S[dim][dim];
    computeStress<dim>(CIJ, ux, S);
    
    //compute residual
    for (unsigned int i=0; i<dim; i++){
      for (unsigned int j=0; j<dim; j++){
	Rux[1,i][1,j] = S[i][j];
      }
    }
    
    //submit residual value for quadrature integration and assemble
    vals.submit_gradient(Rux,q);
  }
}

#endif
