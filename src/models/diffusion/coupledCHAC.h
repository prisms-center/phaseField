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
};

//constructor
template <int dim>
CoupledCHACProblem<dim>::CoupledCHACProblem(): MatrixFreePDE<dim>()
{
  //"n"
  this->getValue["n"]=true; this->getGradient["n"]=true;
  this->setValue["n"]=true; this->setGradient["n"]=true;

  //"c"
  this->getValue["c"]=true; this->getGradient["c"]=true;
  this->setValue["c"]=false; this->setGradient["c"]=true;
}


template <int dim>
void  CoupledCHACProblem<dim>::getRHS(std::map<std::string, typeScalar*>  valsScalar, \
				      std::map<std::string, typeVector*>  valsVector, \
				      unsigned int q) const{
  //parameters
  double Mn[numStructuralOrderParameters]=MnV;
  double Kn[numStructuralOrderParameters]=KnV;
 
 //"n" fields
  scalarvalueType n = valsScalar["n"]->get_value(q);
  scalargradType nx = valsScalar["n"]->get_gradient(q);
  //"c" fields
  scalarvalueType c = valsScalar["c"]->get_value(q);
  scalargradType cx = valsScalar["c"]->get_gradient(q);

  //check to ensure we are working on the intended field
  //order parameter field
  if (this->fields[this->currentFieldIndex].name.compare("n")==0){
    //constants
    scalarType dt=make_vectorized_array(this->timeStep);
    scalarType constN, constNx;
    constN=-Mn[0]; constNx=-Mn[0]*Kn[0]*dt;

    //compute residuals
    valsScalar["n"]->submit_value(constN*rnV,q);
    valsScalar["n"]->submit_gradient(constNx*rnxV,q);
  }
  //concentration field
  else if (this->fields[this->currentFieldIndex].name.compare("c")==0){
    //constants
    scalarType dt=make_vectorized_array(this->timeStep);
    scalarType constCx;
    constCx=-McV*dt;

    //compute residuals
    valsScalar["c"]->submit_gradient(constCx*rcxV,q);
  }
}

#endif
