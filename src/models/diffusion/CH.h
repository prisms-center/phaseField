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
  //RHS implementation for explicit solve
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
  //check if all required parameters correctly specified
#if numFields!=2
#error Compile ERROR: numFields!=2. Number of fields in Cahn-Hilliard problem should be equal to 2.
#endif
#if !defined(McV) || !defined(KcV)
#error Compile ERROR: missing Cahn-Hilliard parameters. Required parameters are McV (mobility) and KcV (length scale parameter).
#endif
#if !defined(rmuV) || !defined(rmuxV) || !defined(rcV) || !defined(rcxV) 
#error Compile ERROR: missing Cahn-Hilliard residual expressions. Required expressions are rmuV, rmuxV, rcV, rcxV
#endif

  //"c"
  this->getValue["c"]=true; this->getGradient["c"]=true;
  this->setValue["c"]=true; this->setGradient["c"]=true;
  //"mu"
  this->getValue["mu"]=true; this->getGradient["mu"]=true;
  this->setValue["mu"]=true; this->setGradient["mu"]=true;
}


template <int dim>
void  CahnHilliardProblem<dim>::getRHS(std::map<std::string, typeScalar*>  valsScalar, \
				       std::map<std::string, typeVector*>  valsVector, \
				       unsigned int q) const{
  //"mu"fields
  scalarvalueType mu = valsScalar["mu"]->get_value(q);
  scalargradType mux = valsScalar["mu"]->get_gradient(q);

  //"c" fields
  scalarvalueType c = valsScalar["c"]->get_value(q);
  scalargradType cx = valsScalar["c"]->get_gradient(q);
  
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
}

#endif
