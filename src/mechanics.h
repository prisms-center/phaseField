//Matrix Free implementation of infinitesimal strain mechanics
#ifndef MECHANICS_H
#define MECHANICS_H

//general headers
#include <fstream>
#include <sstream>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/identity_matrix.h>

//Code specific initializations.  
typedef dealii::parallel::distributed::Vector<double> vectorType;
#if problemDIM==1
typedef dealii::VectorizedArray<double> gradType;
#else 
typedef dealii::Tensor<1, problemDIM, dealii::VectorizedArray<double> > gradType;
#endif

//material models
#include "elasticityModels.h"
#include "computeStress.h"
 

using namespace dealii;

//
//Mechanics class
//
template <int dim>
class MechanicsProblem:  public Subscriptor
{
public:
  MechanicsProblem ();
  void run ();
  void vmult (vectorType &dst, const vectorType &src) const;
private:
  void setup_system ();
  void solve ();
  void output_results (const unsigned int cycle);
  parallel::distributed::Triangulation<dim>      triangulation;
  FESystem<dim>                    fe;
  DoFHandler<dim>                  dof_handler;
  ConstraintMatrix                 constraints;
  IndexSet                         locally_relevant_dofs;
  vectorType                       U, residualU;
  vectorType                       R, H, P;
  double                           setup_time;
  unsigned int                     increment;
  ConditionalOStream               pcout;

  //elasticity matrix
  Table<2, double> CIJ;

  //matrix free objects
  MatrixFree<dim,double>      data;
  vectorType invM;
  void mark_boundaries();
  void  apply_dirichlet_bcs();
  void setup_matrixfree ();
  void updateRHS();
  void computeRHS(const MatrixFree<dim,double> &data, 
		  vectorType &dst, 
		  const vectorType &src,
		  const std::pair<unsigned int,unsigned int> &cell_range) const;
  bool updateRHSValue;
};


//constructor
template <int dim>
MechanicsProblem<dim>::MechanicsProblem ()
  :
  Subscriptor(),
  triangulation (MPI_COMM_WORLD),
  fe (FE_Q<dim>(QGaussLobatto<1>(finiteElementDegree+1)),dim),
  dof_handler (triangulation),
  pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
  CIJ(2*dim-1+dim/3,2*dim-1+dim/3)
{
  double materialConstants[]=MaterialConstantsv;
  getCIJMatrix<dim>(MaterialModelv, materialConstants, CIJ, pcout);
}

//RHS and Matrix-vector product computation
template <int dim>
void MechanicsProblem<dim>::computeRHS (const MatrixFree<dim,double>  &data,
					  vectorType &dst, 
					  const vectorType &src,
					  const std::pair<unsigned int,unsigned int> &cell_range) const
{
  //initialize vals vectors
  FEEvaluation<dim,finiteElementDegree,finiteElementDegree+1,dim,double> vals(data);
  
  //loop over all "cells"
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    
    //read U values for this cell
    vals.reinit (cell); 
    if (updateRHSValue){
      vals.read_dof_values_plain(src);  
    }
    else{
      vals.read_dof_values(src); 
    }
    vals.evaluate (false,true,false); 
    
    //loop over quadrature points
    for (unsigned int q=0; q<vals.n_q_points; ++q){
      Tensor<1, dim, gradType> ux = vals.get_gradient(q);
      Tensor<1, dim, gradType> Cux;
      //compute stress
      dealii::VectorizedArray<double> R[dim][dim];
      computeStress<dim>(CIJ, ux, R);

      //fill residual
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  if (updateRHSValue) {	
	    Cux[1,i][1,j] = -R[i][j];
	  }
	  else{
	    Cux[1,i][1,j] =  R[i][j];
	  }
	}
      }
      //compute residuals
      vals.submit_gradient(Cux,q);
    }
    vals.integrate(false, true); 
    vals.distribute_local_to_global(dst);
  }
}

template <int dim>
void MechanicsProblem<dim>::updateRHS (){
  updateRHSValue=true; 
  //initialize residuals to zero
  residualU=0.0;
  //loop over all cells to compute residuals
  data.cell_loop (&MechanicsProblem::computeRHS, this, residualU, U);
  updateRHSValue=false;
}

// Matrix free data structure vmult operations.
template <int dim>
void MechanicsProblem<dim>::vmult (vectorType &dst, const vectorType &src) const{
  dst=0.0;
  data.cell_loop (&MechanicsProblem::computeRHS, this, dst, src);

  //Account for dirichlet BC's
  const std::vector<unsigned int> &
    constrained_dofs = data.get_constrained_dofs();
  for (unsigned int i=0; i<constrained_dofs.size(); ++i)
    dst(constrained_dofs[i]) += src(constrained_dofs[i]);
}

//setup matrixfree data structures
template <int dim>
void MechanicsProblem<dim>::setup_matrixfree (){
  //setup the matrix free object
  typename MatrixFree<dim,double>::AdditionalData additional_data;
  additional_data.mpi_communicator = MPI_COMM_WORLD;
  additional_data.tasks_parallel_scheme = MatrixFree<dim,double>::AdditionalData::partition_partition;
  additional_data.mapping_update_flags = (update_gradients | update_JxW_values);
  QGaussLobatto<1> quadrature (finiteElementDegree+1);
  data.reinit (dof_handler, constraints, quadrature, additional_data);

  //compute  invM
  data.initialize_dof_vector (invM);
  gradType one;
#if problemDIM==1
  one = make_vectorized_array (1.0);
#else
  for (unsigned int i=0; i<dim; i++)
    one[1,i]=make_vectorized_array (1.0);
#endif

  //select gauss lobatto quad points which are suboptimal but give diogonal M 
  FEEvaluationGL<dim,finiteElementDegree,dim,double> fe_eval(data);
  const unsigned int            n_q_points = fe_eval.n_q_points;
  for (unsigned int cell=0; cell<data.n_macro_cells(); ++cell)
    {
      fe_eval.reinit(cell);
      for (unsigned int q=0; q<n_q_points; ++q)
	fe_eval.submit_value(one,q);
      fe_eval.integrate (true,false);
      fe_eval.distribute_local_to_global (invM);
    }
  invM.compress(VectorOperation::add);
  
  //invert mass matrix diagonal elements
  for (unsigned int k=0; k<invM.local_size(); ++k)
    if (std::abs(invM.local_element(k))>1.0e-15)
      invM.local_element(k) = 1./invM.local_element(k);
    else
      invM.local_element(k) = 0;
}

//setup
template <int dim>
void MechanicsProblem<dim>::setup_system ()
{
  Timer time;
  time.start ();
  setup_time = 0;
  //system_matrix.clear();
  dof_handler.distribute_dofs (fe);
  pcout << "Number of global active cells: " << triangulation.n_global_active_cells() << std::endl;
  pcout << "Number of degrees of freedom:  " << dof_handler.n_dofs() << std::endl;
  DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
  constraints.clear();
  constraints.reinit (locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  //apply Dirichlet BC's
  mark_boundaries();
  apply_dirichlet_bcs();
  constraints.close();
  pcout << "num of constraints:" << constraints.n_constraints() << "\n";
  setup_time += time.wall_time();
  pcout << "Distribute DoFs & B.C.      "
	<<  time.wall_time() << "s" << std::endl;
  time.restart();

  //initialize matrix free data structure
  setup_matrixfree();

  //data structures
  data.initialize_dof_vector(U);
  residualU.reinit(U); 
  R.reinit(U); 
  P.reinit(U);
  H.reinit(U);  

  //initial Condition
  U=0.0;
  
  //timing
  setup_time += time.wall_time();
  pcout << "Wall time for matrix-free setup:    "
	<< time.wall_time() << "s" << std::endl;
}

//solve
template <int dim>
void MechanicsProblem<dim>::solve ()
{
  Timer time; 
  updateRHS(); 
  //cgSolve(U,residualU); 
  SolverControl solver_control(2000, 1e-10*residualU.l2_norm());
  SolverCG<vectorType> cg(solver_control);
  try{
    cg.solve(*this, U, residualU, IdentityMatrix(U.size()));
  }
  catch (...) {
    pcout << "\nWarning: solver did not converge as per set tolerances. consider increasing nsteps or changing tolerances.\n";
  }
  char buffer[200];
  sprintf(buffer, "initial residual:%12.6e, current residual:%12.6e, nsteps:%u, tolerance criterion:%12.6e\n", solver_control.initial_value(), solver_control.last_value(), solver_control.last_step(), solver_control.tolerance()); 
  pcout<<buffer; 
  pcout << "solve time: " << time.wall_time() << "s\n";
}
  
//output 
template <int dim>
void MechanicsProblem<dim>::output_results (const unsigned int cycle)
{
  constraints.distribute (U);
  U.update_ghost_values();
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  //data_out.add_data_vector (U, "u");
  std::vector<std::string> solution_names (dim, "U");
  std::vector<DataComponentInterpretation::DataComponentInterpretation> dataType(dim, DataComponentInterpretation::component_is_part_of_vector);
  data_out.add_data_vector (U, solution_names,
                            DataOut<dim>::type_dof_data,
                            dataType); 
  data_out.build_patches ();
  
  //write to results file
  const std::string filename = "solution-" + Utilities::int_to_string (cycle, std::ceil(std::log10(numIncrements))+1);
  std::ofstream output ((filename +
                         "." + Utilities::int_to_string (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD),std::ceil(std::log10(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)))+1) + ".vtu").c_str());
  data_out.write_vtu (output);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      std::vector<std::string> filenames;
      for (unsigned int i=0;i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i)
	filenames.push_back ("solution-" +
			     Utilities::int_to_string (cycle, std::ceil(std::log10(numIncrements))+1) + "." +
			     Utilities::int_to_string (i, std::ceil(std::log10(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)))+1) + ".vtu");
      std::ofstream master_output ((filename + ".pvtu").c_str());
      data_out.write_pvtu_record (master_output, filenames);
    }
  pcout << "Output written to:" << filename.c_str() << "\n\n" << std::flush;
}


//run
template <int dim>
void MechanicsProblem<dim>::run ()
{
  Timer time;
#if problemDIM==3
  GridGenerator::hyper_rectangle (triangulation, Point<dim>(), Point<dim>(spanX,spanY,spanZ));
#elif problemDIM==2
  GridGenerator::hyper_rectangle (triangulation, Point<dim>(), Point<dim>(spanX,spanY));
#elif problemDIM==1
  GridGenerator::hyper_rectangle (triangulation, Point<dim>(), Point<dim>(spanX));
#endif
  triangulation.refine_global (refineFactor);
  setup_system();

  //write initial conditions to file
  if (writeOutput) output_results(0);

  //Loop over time steps
  for (increment=1; increment<=numIncrements; ++increment)
    {
      pcout << "\nTime increment:" << increment  << ", Wall time: " << time.wall_time() << "s\n" << std::flush;
      //call solve method
      solve ();
      //write results to file
      if ((writeOutput) && (increment%skipOutputSteps==0)) output_results (increment);
    };
  pcout << "Total  time    "  << time.wall_time() << "s" << std::endl; 
  pcout << "PerInc time    "  << time.wall_time()/numIncrements << "s" << std::endl; 
}

#endif
