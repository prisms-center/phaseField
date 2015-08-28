//Matrix Free implementation of the Fickian diffusion model
#ifndef FICKIAN_DIFFUSION_H
#define FICKIAN_DIFFUSION_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

template <int dim>
class FickianProblem: public MatrixFreePDE<dim>
{
 public: 
  FickianProblem();

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
FickianProblem<dim>::FickianProblem(): MatrixFreePDE<dim>()
{
  //check if all required parameters correctly specified
#if numFields!=1
#error Compile ERROR: numFields!=1. Number of fields in Fickian diffusion problem should be equal to 1.
#endif
#if !defined(DcV)
#error Compile ERROR: missing Fickian diffusion parameters. Required parameters are DcV (diffusion coefficient)
#endif
#if !defined(rcV) || !defined(rcxV) 
#error Compile ERROR: missing Fickian diffusion residual expressions. Required expressions are rcV, rcxV
#endif

  //"c"
  this->getValue["c"]=true; this->getGradient["c"]=true;
  this->setValue["c"]=true; this->setGradient["c"]=true;
}


template <int dim>
void  FickianProblem<dim>::getRHS(std::map<std::string, typeScalar*>  valsScalar, \
				       std::map<std::string, typeVector*>  valsVector, \
				       unsigned int q) const{
  //geometric data about the quadrature point
  dealii::Point<dim, dealii::VectorizedArray<double> > point=valsScalar["c"]->quadrature_point(q);
  double x=point[0][0], y=point[1][0];
  
  //time data
  double t=this->currentTime;

  //"c" fields
  scalarvalueType c = valsScalar["c"]->get_value(q);
  scalargradType cx = valsScalar["c"]->get_gradient(q);
  
  //check to ensure we are working on the intended field
  //concentration field
  if (this->fields[this->currentFieldIndex].name.compare("c")==0){
    //compute residuals
    valsScalar["c"]->submit_value(rcV,q);   
    valsScalar["c"]->submit_gradient(rcxV,q);
  }
}

#endif
