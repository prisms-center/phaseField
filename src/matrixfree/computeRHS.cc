//computeRHS() method for MatrixFreePDE class

#ifndef COMPUTERHS_MATRIXFREE_H
#define COMPUTERHS_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//update RHS of each field
template <int dim>
void MatrixFreePDE<dim>::updateRHS(){
  //clear residual vectors before update
  for (unsigned int i=0; i<residualSet.size(); i++){
    (*residualSet[i])=0.0;
  }
  //assembly 
  matrixFreeObject.cell_loop (&MatrixFreePDE<dim>::computeRHS, this, residualSet, solutionSet);
}

//compute RHS
template <int dim>
void  MatrixFreePDE<dim>::computeRHS(const MatrixFree<dim,double> &data, 
				     std::vector<vectorType*> &dst, 
				     const std::vector<vectorType*> &src,
				     const std::pair<unsigned int,unsigned int> &cell_range) const{
  //initialize vals vectors
  std::map<std::string, typeScalar*>  valsScalar;
  std::map<std::string, typeVector*>  valsVector;
  std::map<std::string, unsigned int> valsIndex;
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    if (fields[fieldIndex].type==SCALAR){
      valsScalar[fields[fieldIndex].name]=new typeScalar(data,fieldIndex);
    }
    else{
      valsVector[fields[fieldIndex].name]=new typeVector(data,fieldIndex);
    }
    valsIndex[fields[fieldIndex].name]=fieldIndex;
  }  
  
  //loop over all "cells"
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    unsigned int n_q_points=0;
    //read values from corresponding solution vectors
    for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
      for (std::map<std::string, typeScalar*>::iterator it=valsScalar.begin(); it!=valsScalar.end(); ++it){
	it->second->reinit(cell);
	it->second->read_dof_values_plain(*src[valsIndex[it->first]]);
	it->second->evaluate( getValue.find(it->first)->second,		\
			      getGradient.find(it->first)->second,	\
			      false);
	n_q_points=it->second->n_q_points;
      }
      for (std::map<std::string, typeVector*>::iterator it=valsVector.begin(); it!=valsVector.end(); ++it){
	it->second->reinit(cell);
	it->second->read_dof_values_plain(*src[valsIndex[it->first]]); 
	it->second->evaluate( getValue.find(it->first)->second,		\
			      getGradient.find(it->first)->second,	\
			      false);
	n_q_points=it->second->n_q_points;
      }
    }
    
    //loop over quadrature points
    for (unsigned int q=0; q<n_q_points; ++q){
      getRHS(valsScalar, valsVector, q);
    }

    //Integrate and assemble
    for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
      for (std::map<std::string, typeScalar*>::iterator it=valsScalar.begin(); it!=valsScalar.end(); ++it){
	it->second->integrate(setValue.find(it->first)->second,		\
			      setGradient.find(it->first)->second); 
	it->second->distribute_local_to_global(*dst[valsIndex[it->first]]); 
      }
      for (std::map<std::string, typeVector*>::iterator it=valsVector.begin(); it!=valsVector.end(); ++it){
	it->second->integrate(setValue.find(it->first)->second,		\
			      setGradient.find(it->first)->second); 
	it->second->distribute_local_to_global(*dst[valsIndex[it->first]]); 
      }
    }
  }
}

#endif

