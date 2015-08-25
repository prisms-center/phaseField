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
    //constants
    scalarType dt=make_vectorized_array(this->timeStep);
    scalarType constNx;
    constNx=-KnV*MnV*dt; 
    
    //compute residuals
    valsScalar["n"]->submit_value(rnV,q); 
    valsScalar["n"]->submit_gradient(constNx*rnxV,q);
  }
}

#endif
