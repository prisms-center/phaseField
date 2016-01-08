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
  void getRHS(const MatrixFree<dim,double> &data, 
	      std::vector<vectorType*> &dst, 
	      const std::vector<vectorType*> &src,
	      const std::pair<unsigned int,unsigned int> &cell_range) const;

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
}


template <int dim>
void  FickianProblem<dim>::getRHS(const MatrixFree<dim,double> &data, 
				  std::vector<vectorType*> &dst, 
				  const std::vector<vectorType*> &src,
				  const std::pair<unsigned int,unsigned int> &cell_range) const{
  
  //initialize fields
  typeScalar cVals(data, 0);

  //time data
  double t=this->currentTime;

  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //initialize c field
    cVals.reinit(cell); cVals.read_dof_values_plain(*src[0]); cVals.evaluate(true, true, false);

    //loop over quadrature points
    for (unsigned int q=0; q<cVals.n_q_points; ++q){
      //geometric data about the quadrature point
      dealii::Point<dim, dealii::VectorizedArray<double> > point=cVals.quadrature_point(q);
      double x=point[0][0], y=point[1][0];
      
      //c
      scalarvalueType c = cVals.get_value(q);
      scalargradType cx = cVals.get_gradient(q);
    
      //submit residual value
      cVals.submit_value(rcV,q);
      cVals.submit_gradient(rcxV,q);
    }
    
    //integrate values
    cVals.integrate(true, true);  cVals.distribute_local_to_global(*dst[0]);
  }
}

#endif
