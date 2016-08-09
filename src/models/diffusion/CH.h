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
  void getRHS(const MatrixFree<dim,double> &data, 
	      std::vector<vectorType*> &dst, 
	      const std::vector<vectorType*> &src,
	      const std::pair<unsigned int,unsigned int> &cell_range) const;

  //method to apply initial conditions
  void applyInitialConditions();
  
  //AMR method
  void adaptiveRefine(unsigned int currentIncrement);
  void adaptiveRefineCriterion();
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
}

template <int dim>
void CahnHilliardProblem<dim>::getRHS(const MatrixFree<dim,double> &data, 
				      std::vector<vectorType*> &dst, 
				      const std::vector<vectorType*> &src,
				      const std::pair<unsigned int,unsigned int> &cell_range) const{
  
  //initialize fields
  typeScalar muVals(data, 0), cVals(data,1);

  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //initialize mu field
    muVals.reinit(cell); muVals.read_dof_values_plain(*src[0]); muVals.evaluate(true, true, false);
    
    //initialize c field
    cVals.reinit(cell); cVals.read_dof_values_plain(*src[1]); cVals.evaluate(true, true, false);
    
    //loop over quadrature points
    for (unsigned int q=0; q<cVals.n_q_points; ++q){
      //mu
      scalarvalueType mu = muVals.get_value(q);
      scalargradType mux = muVals.get_gradient(q);
      
      //c
      scalarvalueType c = cVals.get_value(q);
      scalargradType cx = cVals.get_gradient(q);
      
      //submit values
      muVals.submit_value(rmuV,q); muVals.submit_gradient(rmuxV,q);
      cVals.submit_value(rcV,q); cVals.submit_gradient(rcxV,q);
    }
    
    //integrate values
    muVals.integrate(true, true);  muVals.distribute_local_to_global(*dst[0]);
    cVals.integrate(true, true);  cVals.distribute_local_to_global(*dst[1]);
  }
}

#ifndef hAdaptivity
//adaptive refinement control
template <int dim>
void CahnHilliardProblem<dim>::adaptiveRefine(unsigned int currentIncrement){
}

//adaptive refinement criterion
template <int dim>
void CahnHilliardProblem<dim>::adaptiveRefineCriterion(){
}
#endif

#endif
