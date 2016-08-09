//Matrix Free implementation of infinitesimal strain mechanics
#ifndef MECHANICS_MECHANICS_H
#define MECHANICS_MECHANICS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

//material models
#include "computeStress.h"

template <int dim>
class MechanicsProblem: public MatrixFreePDE<dim>
{
 public: 
  MechanicsProblem();

 private:
  //elasticity matrix
  Table<2, double> CIJ;
  
  //RHS implementation for explicit solve
  void getRHS(const MatrixFree<dim,double> &data, 
	      std::vector<vectorType*> &dst, 
	      const std::vector<vectorType*> &src,
	      const std::pair<unsigned int,unsigned int> &cell_range) const;
    
  //LHS implementation for implicit solve 
  void  getLHS(const MatrixFree<dim,double> &data, 
	       vectorType &dst, 
	       const vectorType &src,
	       const std::pair<unsigned int,unsigned int> &cell_range) const;
  
  //methods to apply dirichlet BC's
  void markBoundaries();
  void applyDirichletBCs();
  
  //AMR method
  void adaptiveRefine(unsigned int currentIncrement);
  void adaptiveRefineCriterion();
}; 

//constructor
template <int dim>
MechanicsProblem<dim>::MechanicsProblem(): MatrixFreePDE<dim>(), 
  CIJ(2*dim-1+dim/3,2*dim-1+dim/3)
{
  //check if all required parameters correctly specified
#if numFields!=problemDIM
#error Compile ERROR: numFields!=problemDIM. Number of fields in Mechanics problem should be equal to the problem dimension
#endif
#ifndef MaterialModelV
#error Compile ERROR: missing material property variable: MaterialModelV
#endif
#ifndef MaterialConstantsV
#error Compile ERROR: missing material property variable: MaterialConstantsV
#endif

  //initialize elasticity matrix
#if defined(MaterialModelV) && defined(MaterialConstantsV)
  double materialConstants[]=MaterialConstantsV;
  getCIJMatrix<dim>(MaterialModelV, materialConstants, CIJ, this->pcout);
#endif
}

//implementation of the RHS evaluation 
template <int dim>
void  MechanicsProblem<dim>::getRHS(const MatrixFree<dim,double> &data, 
				    std::vector<vectorType*> &dst, 
				    const std::vector<vectorType*> &src,
				    const std::pair<unsigned int,unsigned int> &cell_range) const{
 
  typeVector uVals(data, 0);
  
  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //initialize u field 
    uVals.reinit(cell); uVals.read_dof_values_plain(*src[0]); uVals.evaluate(false, true, false);

    //loop over quadrature points
    for (unsigned int q=0; q<uVals.n_q_points; ++q){
      //u
      vectorgradType ux = uVals.get_gradient(q);
      vectorgradType Rux;

      //compute strain tensor
      dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]);
	}
      }
      
      //compute stress tensor
      computeStress<dim>(CIJ, E, S);

      //compute residual
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  Rux[i][j] -= S[i][j]; 
	}
      }
      
      //submit residual value
      uVals.submit_gradient(Rux,q);
    }
    
    //integrate
    uVals.integrate(false, true); uVals.distribute_local_to_global(*dst[0]);
  }
}

template <int dim>
void  MechanicsProblem<dim>::getLHS(const MatrixFree<dim,double> &data, 
				    vectorType &dst, 
				    const vectorType &src,
				    const std::pair<unsigned int,unsigned int> &cell_range) const{
  typeVector uVals(data, 0);
  
  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //initialize u field 
    uVals.reinit(cell); uVals.read_dof_values_plain(src); uVals.evaluate(false, true, false);

    //loop over quadrature points
    for (unsigned int q=0; q<uVals.n_q_points; ++q){
      //u
      vectorgradType ux = uVals.get_gradient(q);
      vectorgradType Rux;

      //compute strain tensor
      dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]);
	}
      }
    
      //compute stress tensor
      computeStress<dim>(CIJ, E, S);
      
      //compute residual
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  Rux[i][j] = S[i][j]; 
	}
      }
      
      //submit residual value
      uVals.submit_gradient(Rux,q);
    }
    
    //integrate
    uVals.integrate(false, true); uVals.distribute_local_to_global(dst);
  }
}

#ifndef hAdaptivity
//adaptive refinement control
template <int dim>
void MechanicsProblem<dim>::adaptiveRefine(unsigned int currentIncrement){
}

//adaptive refinement criterion
template <int dim>
void MechanicsProblem<dim>::adaptiveRefineCriterion(){
}
#endif

#endif

