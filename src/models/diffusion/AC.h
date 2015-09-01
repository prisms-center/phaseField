//Matrix Free implementation of Allen-Cahn order parameter dynamics
#ifndef AC_DIFFUSION_H
#define AC_DIFFUSION_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

template <int dim>
class AllenCahnProblem: public MatrixFreePDE<dim>
{
 public: 
  AllenCahnProblem();

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
AllenCahnProblem<dim>::AllenCahnProblem(): MatrixFreePDE<dim>()
{
  //check if all required parameters correctly specified
#if numFields!=1
#error Compile ERROR: numFields!=1. Number of fields in standard Allen-Cahn problem should be equal to 1.
#endif
#if !defined(MnV) || !defined(KnV)
#error Compile ERROR: missing Allen-Cahn parameters. Required parameters are MnV (mobility) and KnV (length scale parameter).
#endif
#if !defined(rnV) || !defined(rnxV)
#error Compile ERROR: missing Allen-Cahn residual expressions. Required expressions are rnV, rnxV.
#endif

  //"n"
  this->getValue["n"]=true; this->getGradient["n"]=true;
  this->setValue["n"]=true; this->setGradient["n"]=true;
}


template <int dim>
void  AllenCahnProblem<dim>::getRHS(std::map<std::string, typeScalar*>  valsScalar, \
				    std::map<std::string, typeVector*>  valsVector, \
				    unsigned int q) const{
  //"n"fields
  scalarvalueType n = valsScalar["n"]->get_value(q);
  scalargradType nx = valsScalar["n"]->get_gradient(q);
  
  //check to ensure we are working on the intended field
  //order parameter field
  if (this->fields[this->currentFieldIndex].name.compare("n")==0){
    //compute residuals
    valsScalar["n"]->submit_value(rnV,q); 
    valsScalar["n"]->submit_gradient(rnxV,q);
  }
}

#endif
