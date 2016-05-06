//Matrix Free implementation of coupled Cahn-Hilliard and Mechanics formulation 
#ifndef CHMECHANICS_DIFFUSION_H
#define CHMECHANICS_DIFFUSION_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#include "../../../include/matrixFreePDE.h"

//material models
#include "../mechanics/computeStress.h"

template <int dim>
class CoupledCahnHilliardMechanicsProblem: public MatrixFreePDE<dim>
{
 public: 
  CoupledCahnHilliardMechanicsProblem();

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
  
  //method to apply initial conditions
  void applyInitialConditions();
 
  //methods to apply dirichlet BC's on displacement
  void applyDirichletBCs();


  void getEnergy(const MatrixFree<dim,double> &data,
  				    std::vector<vectorType*> &dst,
  				    const std::vector<vectorType*> &src,
  				    const std::pair<unsigned int,unsigned int> &cell_range);
  Threads::Mutex assembler_lock;
};

//constructor
template <int dim>
CoupledCahnHilliardMechanicsProblem<dim>::CoupledCahnHilliardMechanicsProblem(): MatrixFreePDE<dim>(),
  CIJ(2*dim-1+dim/3,2*dim-1+dim/3)
{
  //check if all required parameters correctly specified
#if numFields!=(problemDIM+2)
#error Compile ERROR: numFields!=(problemDIM+2). Number of fields in coupled Cahn-Hilliard and Mechanics problem should be equal to (problemDIM+2).
#endif
#if !defined(McV) || !defined(KcV)
#error Compile ERROR: missing Cahn-Hilliard parameters. Required parameters are McV (mobility) and KcV (length scale parameter).
#endif
#if !defined(rmuV) || !defined(rmuxV) || !defined(rcV) || !defined(rcxV) 
#error Compile ERROR: missing Cahn-Hilliard residual expressions. Required expressions are rmuV, rmuxV, rcV, rcxV.
#endif
#ifndef MaterialModelV
#error Compile ERROR: missing material property variable: MaterialModelV.
#endif
#ifndef MaterialConstantsV
#error  Compile ERROR: missing material property variable: MaterialConstantsV.
#endif
#ifndef chemicalStrainV
#error  Compile ERROR: missing coupled property variable: chemicalStrainV.
#endif

  //initialize elasticity matrix
#if defined(MaterialModelV) && defined(MaterialConstantsV)
  double materialConstants[]=MaterialConstantsV;
  getCIJMatrix<dim>(MaterialModelV, materialConstants, CIJ, this->pcout);
#endif
}


template <int dim>
void  CoupledCahnHilliardMechanicsProblem<dim>::getRHS(const MatrixFree<dim,double> &data, 
						       std::vector<vectorType*> &dst, 
						       const std::vector<vectorType*> &src,
						       const std::pair<unsigned int,unsigned int> &cell_range) const{

  //initialize fields
  typeScalar muVals(data, 0), cVals(data,1);
  typeVector uVals(data, 2);

  //define identity tensor
  vectorgradType I; 
  for (unsigned int i=0; i<dim; i++) {I[i][i]+=constV(1.0);}
  
  //loop over cells
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //initialize mu field
    muVals.reinit(cell); muVals.read_dof_values_plain(*src[0]); muVals.evaluate(true, true, false);

    //initialize c field
    cVals.reinit(cell); cVals.read_dof_values_plain(*src[1]); cVals.evaluate(true, true, false);

    //initialize u field 
    uVals.reinit(cell); uVals.read_dof_values_plain(*src[2]); uVals.evaluate(false, true, false);

    //loop over quadrature points
    for (unsigned int q=0; q<cVals.n_q_points; ++q){
      //mu
      scalarvalueType mu = muVals.get_value(q);
      scalargradType mux = muVals.get_gradient(q);
      //c
      scalarvalueType c = cVals.get_value(q);
      scalargradType cx = cVals.get_gradient(q);
      //u
      vectorgradType ux = uVals.get_gradient(q);
      vectorgradType Rux;
      
      //apply chemical strain
      ux+=chemicalStrainV*I;

      //compute strain tensor
      dealii::VectorizedArray<double> E[dim][dim], S[dim][dim];
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  E[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]);
	}
      }
      
      //compute stress tensor
      computeStress<dim>(CIJ, E, S);

      //fill residual corresponding to mechanics
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  Rux[i][j] -= S[i][j];
	}
      }
      
      //submit values
      muVals.submit_value(rmuV,q); muVals.submit_gradient(rmuxV,q);
      cVals.submit_value(rcV,q); cVals.submit_gradient(rcxV,q);
      uVals.submit_gradient(Rux,q);
    }
    
    //integrate values
    muVals.integrate(true, true);  muVals.distribute_local_to_global(*dst[0]);
    cVals.integrate(true, true);   cVals.distribute_local_to_global(*dst[1]);
    uVals.integrate(false, true); uVals.distribute_local_to_global(*dst[2]);
  }
}

template <int dim>
void  CoupledCahnHilliardMechanicsProblem<dim>::getLHS(const MatrixFree<dim,double> &data, 
						       vectorType &dst, 
						       const vectorType &src,
						       const std::pair<unsigned int,unsigned int> &cell_range) const{
  typeVector uVals(data, 2);
  
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

// Calculate the free energy
template <int dim>
void  CoupledCahnHilliardMechanicsProblem<dim>::getEnergy(const MatrixFree<dim,double> &data,
				    std::vector<vectorType*> &dst,
				    const std::vector<vectorType*> &src,
				    const std::pair<unsigned int,unsigned int> &cell_range) {

	//initialize fields
	  typeScalar muVals(data, 0), cVals(data,1);
	  typeVector uVals(data, 2);

	  //loop over cells
	  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){

	    //initialize c field
	    cVals.reinit(cell); cVals.read_dof_values_plain(*src[1]); cVals.evaluate(true, true, false);

	    //initialize u field
	    uVals.reinit(cell); uVals.read_dof_values_plain(*src[2]); uVals.evaluate(false, true, false);

	    dealii::AlignedVector<dealii::VectorizedArray<double> > JxW(cVals.n_q_points);
	    cVals.fill_JxW_values(JxW);

	    //loop over quadrature points
	    for (unsigned int q=0; q<cVals.n_q_points; ++q){
	      //c
	      scalarvalueType c = cVals.get_value(q);
	      scalargradType cx = cVals.get_gradient(q);
	      //u
	      vectorgradType ux = uVals.get_gradient(q);
	      vectorgradType Rux;

	      scalarvalueType total_energy_density = constV(0.0);

	      total_energy_density = c*c*c*c - constV(2.0)*c*c*c + c*c;

	      for (unsigned int i=1; i<dim; i++){
	    	  total_energy_density += KcV*constV(0.5)*cx[i]*cx[i];
	      }

	      assembler_lock.acquire ();
	      for (unsigned i=0; i<2;i++){
	    	  this->energy+=total_energy_density[i]*JxW[q][i];
	      }
	      assembler_lock.release ();
	    }
	}
}


#endif
