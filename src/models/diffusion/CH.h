//Matrix Free implementation of Cahn-Hilliard spinodal decomposition mixed (split) formulation 
#ifndef CH_DIFFUSION_H
#define CH_DIFFUSION_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

template <int dim>
class CahnHilliardProblem: public MatrixFreePDE<dim>
{
 public: 
  CahnHilliardProblem();

 private:
  //RHS implementation for implicit/explicit solve
  void getRHS(std::map<std::string, typeScalar*>  valsScalar, \
	      std::map<std::string, typeVector*>  valsVector, \
	      unsigned int q) const;
  
  //method to apply initial conditions
  void applyInitialConditions();
};

//constructor
template <int dim>
CahnHilliardProblem<dim>::CahnHilliardProblem(): MatrixFreePDE<dim>()
{
  //"c"
  this->getValue["c"]=true; this->getGradient["c"]=true;
  this->setValue["c"]=true; this->setGradient["c"]=true;
  //"mu"
  this->getValue["mu"]=true; this->getGradient["mu"]=true;
  this->setValue["mu"]=true; this->setGradient["mu"]=true;
}

//implementation of the RHS evaluation 
template <int dim>
void  CahnHilliardProblem<dim>::getRHS(std::map<std::string, typeScalar*>  valsScalar, \
				       std::map<std::string, typeVector*>  valsVector, \
				       unsigned int q) const{
  //"c" fields
  scalarvalueType c = valsScalar["c"]->get_value(q);
  scalargradType cx = valsScalar["c"]->get_gradient(q);
  
  //"mu"fields
  scalarvalueType mu = valsScalar["mu"]->get_value(q);
  scalargradType mux = valsScalar["mu"]->get_gradient(q);

  //constants
  scalarType constMux, constCx;
  scalarType dt=make_vectorized_array(this->timeStep);
  constMux=KcV; constCx=-McV*dt;

  //compute residuals
  valsScalar["c"]->submit_value(rcV,q);   valsScalar["c"]->submit_gradient(constCx*rcxV,q);
  valsScalar["mu"]->submit_value(rmuV,q); valsScalar["mu"]->submit_gradient(constMux*rmuxV,q);
}

#endif
