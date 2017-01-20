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
  void getRHS(const MatrixFree<dim,double> &data, 
	      std::vector<vectorType*> &dst, 
	      const std::vector<vectorType*> &src,
	      const std::pair<unsigned int,unsigned int> &cell_range) const;

  // Method to apply initial conditions
  void applyInitialConditions();

  // Method to make modification to one of the fields (used for nucleation)
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
}


template <int dim>
void  CoupledCHACProblem<dim>::getRHS(const MatrixFree<dim,double> &data, 
				      std::vector<vectorType*> &dst, 
				      const std::vector<vectorType*> &src,
				      const std::pair<unsigned int,unsigned int> &cell_range) const{
  
  //initialize fields
  typeScalar nVals(data, 0), cVals(data,1);

  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //initialize mu field
    nVals.reinit(cell); nVals.read_dof_values_plain(*src[0]); nVals.evaluate(true, true, false);
    
    //initialize c field
    cVals.reinit(cell); cVals.read_dof_values_plain(*src[1]); cVals.evaluate(true, true, false);
    
    //loop over quadrature points
    for (unsigned int q=0; q<cVals.n_q_points; ++q){
      //n
      scalarvalueType n = nVals.get_value(q);
      scalargradType nx = nVals.get_gradient(q);
      
      //c
      scalarvalueType c = cVals.get_value(q);
      scalargradType cx = cVals.get_gradient(q);
      
      //submit values
      nVals.submit_value(rnV,q); nVals.submit_gradient(rnxV,q);
      cVals.submit_value(rcV,q); cVals.submit_gradient(rcxV,q);
    }
    
    //integrate values
    nVals.integrate(true, true);  nVals.distribute_local_to_global(*dst[0]);
    cVals.integrate(true, true);  cVals.distribute_local_to_global(*dst[1]);
  }
}


#endif
