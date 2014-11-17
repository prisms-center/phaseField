//Matrix Free implementation of infinitesimal strain mechanics
#ifndef MECHANICS2_H
#define MECHANICS2_H
#include "matrixFreePDE.h"

//material models
#include "../elasticityModels.h"
#include "../computeStress.h"

typedef  FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,1,double>           typeScalar;
typedef  FEEvaluation<problemDIM,finiteElementDegree,finiteElementDegree+1,problemDIM,double>  typeVector;

template <int dim>
class MechanicsProblem: public MatrixFreePDE<dim>
{
 public: 
  MechanicsProblem();
  void computeRHS (const MatrixFree<dim,double> &data, 
		   std::vector<vectorType*> &dst, 
		   const std::vector<vectorType*> &src,
		   const std::pair<unsigned int,unsigned int> &cell_range) const;
 private:
  //elasticity matrix
  Table<2, double> CIJ;
  void getRHS(std::map<std::string, typeScalar*>  valsScalar, \
	      std::map<std::string, typeVector*>  valsVector, \
	      unsigned int q) const;  
};

//constructor
template <int dim>
MechanicsProblem<dim>::MechanicsProblem(): MatrixFreePDE<dim>(), 
  CIJ(2*dim-1+dim/3,2*dim-1+dim/3)
{
  //initialize elasticity matrix
  double materialConstants[]=MaterialConstantsv;
  getCIJMatrix<dim>(MaterialModelv, materialConstants, CIJ, this->pcout);
}


template <int dim>
void  MechanicsProblem<dim>::getRHS(std::map<std::string, typeScalar*>  valsScalar, \
				    std::map<std::string, typeVector*>  valsVector, \
				    unsigned int q) const{
  Tensor<1, dim, gradType> ux = valsVector["u"]->get_gradient(q);
  Tensor<1, dim, gradType> Cux;
  //compute stress
  dealii::VectorizedArray<double> R[dim][dim];
  computeStress<dim>(CIJ, ux, R);
  
  //fill residual
  for (unsigned int i=0; i<dim; i++){
    for (unsigned int j=0; j<dim; j++){
      Cux[1,i][1,j] = -R[i][j];
    }
  }
  
  //compute residuals
  valsVector["u"]->submit_gradient(Cux,q);
}


//compute RHS
template <int dim>
void  MechanicsProblem<dim>::computeRHS(const MatrixFree<dim,double> &data, 
					std::vector<vectorType*> &dst, 
					const std::vector<vectorType*> &src,
					const std::pair<unsigned int,unsigned int> &cell_range) const{
  //initialize vals vectors
  std::map<std::string, typeScalar*>  valsScalar;
  std::map<std::string, typeVector*>  valsVector;
  std::map<std::string, unsigned int> valsIndex;
  for(unsigned int fieldIndex=0; fieldIndex<this->fields.size(); fieldIndex++){
    if (this->fields[fieldIndex].type==SCALAR){
      valsScalar[this->fields[fieldIndex].name]=new typeScalar(data,fieldIndex);
    }
    else{
      valsVector[this->fields[fieldIndex].name]=new typeVector(data,fieldIndex);
    }
    valsIndex[this->fields[fieldIndex].name]=fieldIndex;
  }  
  
  //loop over all "cells"
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    unsigned int n_q_points=0;
    //read values from corresponding solution vectors
    for(unsigned int fieldIndex=0; fieldIndex<this->fields.size(); fieldIndex++){
      for (std::map<std::string, typeScalar*>::iterator it=valsScalar.begin(); it!=valsScalar.end(); ++it){
	it->second->reinit(cell);
	it->second->read_dof_values_plain(*src[valsIndex[it->first]]); 
	it->second->evaluate(true,true,false);
	n_q_points=it->second->n_q_points;
      }
      for (std::map<std::string, typeVector*>::iterator it=valsVector.begin(); it!=valsVector.end(); ++it){
	it->second->reinit(cell);
	it->second->read_dof_values_plain(*src[valsIndex[it->first]]); 
	it->second->evaluate(true,true,false);
	n_q_points=it->second->n_q_points;
      }
    }
    
    //loop over quadrature points
    for (unsigned int q=0; q<n_q_points; ++q){
      getRHS(valsScalar, valsVector, q);
    }

    //Integrate and assemble
    for(unsigned int fieldIndex=0; fieldIndex<this->fields.size(); fieldIndex++){
      for (std::map<std::string, typeScalar*>::iterator it=valsScalar.begin(); it!=valsScalar.end(); ++it){
	it->second->integrate(false, true); 
	it->second->distribute_local_to_global(*dst[valsIndex[it->first]]); 
      }
      for (std::map<std::string, typeVector*>::iterator it=valsVector.begin(); it!=valsVector.end(); ++it){
	it->second->integrate(false, true); 
	it->second->distribute_local_to_global(*dst[valsIndex[it->first]]); 
      }
    }
  }
}

#endif
