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
  void getRHS(const MatrixFree<dim,double> &data, 
	      std::vector<vectorType*> &dst, 
	      const std::vector<vectorType*> &src,
	      const std::pair<unsigned int,unsigned int> &cell_range) const;

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
}


template <int dim>
void  AllenCahnProblem<dim>::getRHS(const MatrixFree<dim,double> &data, 
				    std::vector<vectorType*> &dst, 
				    const std::vector<vectorType*> &src,
				    const std::pair<unsigned int,unsigned int> &cell_range) const{
  
  //initialize fields
  typeScalar nVals(data, 0);

  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //initialize mu field
    nVals.reinit(cell); nVals.read_dof_values_plain(*src[0]); nVals.evaluate(true, true, false);
    
    //loop over quadrature points
    for (unsigned int q=0; q<nVals.n_q_points; ++q){
      //mu
      scalarvalueType n = nVals.get_value(q);
      scalargradType nx = nVals.get_gradient(q);
      
      //submit values
      nVals.submit_value(rnV,q); nVals.submit_gradient(rnxV,q);
    }
    
    //integrate values
    nVals.integrate(true, true);  nVals.distribute_local_to_global(*dst[0]);
  }
}

#endif
