//Matrix Free implementation of coupled Cahn-Hilliard and Allen-Cahn dynamics
#ifndef CHAC_DIFFUSION_H
#define CHAC_DIFFUSION_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

template <int dim>
class CoupledCHACProblem: public MatrixFreePDE<dim>
{
 public: 
  CoupledCHACProblem();

 private:
  //RHS implementation for explicit solve
  void getRHS(std::map<std::string, typeScalar*>  valsScalar, \
	      std::map<std::string, typeVector*>  valsVector, \
	      unsigned int q) const;

  //method to apply initial conditions
  void applyInitialConditions();
  void modifySolutionFields();
};

//constructor
template <int dim>
CoupledCHACProblem<dim>::CoupledCHACProblem(): MatrixFreePDE<dim>()
{
  //check if all required parameters correctly specified
#if numFields!=2
#error Compile ERROR: numFields!=2. Number of fields in coupled Cahn-Hilliard and Allen-Cahn problem should be equal to 2.
#endif
#if !defined(McV) 
#error Compile ERROR: missing Cahn-Hilliard parameters. Required parameters are McV (mobility).
#endif
#if !defined(MnV) || !defined(KnV)  
#error Compile ERROR: missing Allen-Cahn parameters. Required parameters are MnV (mobility) and KnV (length scale parameter).
#endif
#if !defined(rnV) || !defined(rnxV) || !defined(rcV) || !defined(rcxV) 
#error Compile ERROR: missing required residual expressions. Required expressions are rnV, rnxV, rcV, rcxV
#endif

  //"n"
  this->getValue["n"]=true; this->getGradient["n"]=true;
  this->setValue["n"]=true; this->setGradient["n"]=true;

  //"c"
  this->getValue["c"]=true; this->getGradient["c"]=true;
  this->setValue["c"]=true; this->setGradient["c"]=true;
}


template <int dim>
void  CoupledCHACProblem<dim>::getRHS(std::map<std::string, typeScalar*>  valsScalar, \
				      std::map<std::string, typeVector*>  valsVector, \
				      unsigned int q) const{
 //"n" fields
  scalarvalueType n = valsScalar["n"]->get_value(q);
  scalargradType nx = valsScalar["n"]->get_gradient(q);
  //"c" fields
  scalarvalueType c = valsScalar["c"]->get_value(q);
  scalargradType cx = valsScalar["c"]->get_gradient(q);

  //check to ensure we are working on the intended field
  //order parameter field
  if (this->fields[this->currentFieldIndex].name.compare("n")==0){
    //compute residuals
    valsScalar["n"]->submit_value(rnV,q);
    valsScalar["n"]->submit_gradient(rnxV,q);
  }
  //concentration field
  else if (this->fields[this->currentFieldIndex].name.compare("c")==0){
    //compute residuals
    valsScalar["c"]->submit_value(rcV,q);
    valsScalar["c"]->submit_gradient(rcxV,q);
  }
}

#endif
