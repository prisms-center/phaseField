//base class for matrix Free implementation of PDE's
#ifndef MATRIXFREEPDE_H
#define MATRIXFREEPDE_H

//general headers
#include <fstream>
#include <sstream>

//dealii headers
#include "../../include/dealIIheaders.h"

//PRISMS headers
#include "fields.h"

//define data types  
typedef dealii::parallel::distributed::Vector<double> vectorType;
#if problemDIM==1
typedef dealii::VectorizedArray<double> gradType;
#else 
typedef dealii::Tensor<1, problemDIM, dealii::VectorizedArray<double> > gradType;
#endif

using namespace dealii;

//
//base class for matrix free PDE's
//
template <int dim>
class MatrixFreePDE:public Subscriptor
{
 public:
  MatrixFreePDE(); 
  ~MatrixFreePDE(); 
  void init  ();
  void solve ();
  void vmult (vectorType &dst, const vectorType &src) const;
  std::vector<Field<dim> >                  fields;
 private:
  void solveIncrement ();
  void outputResults  ();
  parallel::distributed::Triangulation<dim> triangulation;
  std::vector<FESystem<dim>*>          FESet;
  std::vector<const ConstraintMatrix*> constraintsSet;
  std::vector<const DoFHandler<dim>*>  dofHandlersSet;
  std::vector<const IndexSet*>         locally_relevant_dofsSet;
  std::vector<vectorType*>             solutionSet, residualSet;

  //matrix free objects
  MatrixFree<dim,double>               matrixFreeObject;
  vectorType                           invM;
  
  //matrix free methods
  void updateRHS();
  void computeRHS(const MatrixFree<dim,double> &data, 
		  std::vector<vectorType*> &dst, 
		  const std::vector<vectorType*> &src,
		  const std::pair<unsigned int,unsigned int> &cell_range) const;
  void computeInvM();
  
  //variables for time dependent problems 
  //isTimeDependentBVP flag is used to see if invM, time steppping in
  //run(), etc are necessary
  bool isTimeDependentBVP;
  double timeStep, currentTime, finalTime;
  unsigned int currentIncrement, totalIncrements;
  
  //parallel message stream
  ConditionalOStream  pcout;  
  //compute time log
  TimerOutput computing_timer;
};

 //constructor
 template <int dim>
 MatrixFreePDE<dim>::MatrixFreePDE ()
 :
 Subscriptor(),
   triangulation (MPI_COMM_WORLD),
   isTimeDependentBVP(false),
   currentTime(0.0),
   currentIncrement(0),
   pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
   computing_timer (pcout, TimerOutput::summary, TimerOutput::wall_times)
 {}

 //destructor
 template <int dim>
 MatrixFreePDE<dim>::~MatrixFreePDE ()
 {
   matrixFreeObject.clear();
   for(unsigned int iter=0; iter<fields.size(); iter++){
     delete locally_relevant_dofsSet[iter];
     delete constraintsSet[iter];
     delete dofHandlersSet[iter];
     delete FESet[iter];
     delete solutionSet[iter];
     delete residualSet[iter];
   } 
 }

 //populate with fields and setup matrix free system
 template <int dim>
 void MatrixFreePDE<dim>::init(){
 //void MatrixFreePDE<dim>::init(std::vector<Field<dim> >& _fields){
   computing_timer.enter_section("matrixFreePDE: initialization"); 

   //creating mesh
   pcout << "creating problem mesh\n";
 #if problemDIM==3
   GridGenerator::hyper_rectangle (triangulation, Point<dim>(), Point<dim>(spanX,spanY,spanZ));
 #elif problemDIM==2
   GridGenerator::hyper_rectangle (triangulation, Point<dim>(), Point<dim>(spanX,spanY));
 #elif problemDIM==1
   GridGenerator::hyper_rectangle (triangulation, Point<dim>(), Point<dim>(spanX));
 #endif
   triangulation.refine_global (refineFactor);
   pcout << "number of elements: " << triangulation.n_global_active_cells() << std::endl;

   //setup system
   unsigned int totalDOFs=0;
   for(typename std::vector<Field<dim> >::iterator it = fields.begin(); it != fields.end(); ++it){
     //print to std::out
     char buffer[100];
     sprintf(buffer,"initializing finite element space P^%u for %6s:%9s field '%s'\n", \
	     finiteElementDegree,					\
	     (it->type==SCALAR ? "SCALAR":"VECTOR"),			\
	     (it->pdetype==PARABOLIC ? "PARABOLIC":"ELLIPTIC"), it->name.c_str());
     pcout << buffer;

     //check if any time dependent fields present
     if (it->pdetype==PARABOLIC){
       isTimeDependentBVP=true;
     }

     //create FESystem
     FESystem<dim>* fe;
     if (it->type==SCALAR){
       fe=new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(finiteElementDegree+1)),1);
     }
     else if (it->type==VECTOR){
       fe=new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(finiteElementDegree+1)),dim);
     }
     else{
       pcout << "\nmatrixFreePDE.h: unknown field type\n";
       exit(-1);
     }
     FESet.push_back(fe);

     //distribute DOFs
     DoFHandler<dim>* dof_handler=new DoFHandler<dim>(triangulation);
     dof_handler->distribute_dofs (*fe);
     dofHandlersSet.push_back(dof_handler); 
     totalDOFs+=dof_handler->n_dofs();

     //extract locally_relevant_dofs
     IndexSet* locally_relevant_dofs=new IndexSet;
     DoFTools::extract_locally_relevant_dofs (*dof_handler, *locally_relevant_dofs);
     locally_relevant_dofsSet.push_back(locally_relevant_dofs);

     //create constraints
     ConstraintMatrix* constraints=new ConstraintMatrix;
     constraints->clear();
     constraints->reinit(*locally_relevant_dofs);
     DoFTools::make_hanging_node_constraints (*dof_handler, *constraints);
     //apply zero Dirichlet BC's for ELLIPTIC fields. This is just the
     //default and can be changed later in the specific BVP
     //implementation
     if (it->pdetype==ELLIPTIC){
       VectorTools::interpolate_boundary_values (*dof_handler, 0, ZeroFunction<dim>(it->numComponents), *constraints);
     }
     constraints->close();  
     constraintsSet.push_back(constraints);
   }
   pcout << "number of degrees of freedom: " << totalDOFs << std::endl;

   //setup the matrix free object
   typename MatrixFree<dim,double>::AdditionalData additional_data;
   additional_data.mpi_communicator = MPI_COMM_WORLD;
   additional_data.tasks_parallel_scheme = MatrixFree<dim,double>::AdditionalData::partition_partition;
   additional_data.mapping_update_flags = (update_values | update_gradients | update_JxW_values);
   QGaussLobatto<1> quadrature (finiteElementDegree+1);
   matrixFreeObject.reinit (dofHandlersSet, constraintsSet, quadrature, additional_data);
   pcout << "completed initialization of the matrix free object\n";
 
   //setup problem vectors
   pcout << "initializing parallel::distributed residual and solution vectors\n";
   for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
     vectorType* U=new vectorType;
     vectorType* R=new vectorType;
     matrixFreeObject.initialize_dof_vector(*U,  fieldIndex);
     matrixFreeObject.initialize_dof_vector(*R,  fieldIndex);
     solutionSet.push_back(U);
     residualSet.push_back(R);
   }

   //check if time dependent BVP and compute invM
   if (isTimeDependentBVP){
     computeInvM();
   }
   computing_timer.exit_section("matrixFreePDE: initialization");  
}

//compute inverse of the diagonal mass matrix and store in vector invM
template <int dim>
void MatrixFreePDE<dim>::computeInvM(){
  //initialize  invM
  bool invMInitialized=false;
  unsigned int parabolicFieldIndex=0;
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    if (fields[fieldIndex].pdetype==PARABOLIC){
      matrixFreeObject.initialize_dof_vector (invM, fieldIndex);
      parabolicFieldIndex=fieldIndex;
      invMInitialized=true;
      break;
    }
  }
  //check if invM initialized
  if (!invMInitialized){
    pcout << "matrixFreePDE.h: no PARABOLIC field... hence cannot initialize invM\n";
    exit(-1);
  }
  
  //compute invM
  VectorizedArray<double> one = make_vectorized_array (1.0);
  
  //select gauss lobatto quad points which are suboptimal but give diogonal M 
  FEEvaluation<dim,finiteElementDegree>  fe_eval(matrixFreeObject, parabolicFieldIndex);
  const unsigned int n_q_points = fe_eval.n_q_points;
  for (unsigned int cell=0; cell<matrixFreeObject.n_macro_cells(); ++cell){
    fe_eval.reinit(cell);
    for (unsigned int q=0; q<n_q_points; ++q){
      fe_eval.submit_value(one,q);
    }
    fe_eval.integrate (true,false);
    fe_eval.distribute_local_to_global (invM);
  }
  invM.compress(VectorOperation::add);
  
  //invert mass matrix diagonal elements
  for (unsigned int k=0; k<invM.local_size(); ++k){
    if (std::abs(invM.local_element(k))>1.0e-15){
      invM.local_element(k) = 1./invM.local_element(k);
    }
    else{
      invM.local_element(k) = 0;
    }
  } 
  pcout << "computed mass matrix (using FE space for field: " << parabolicFieldIndex << ")\n";
}

//compute RHS
template <int dim>
void  MatrixFreePDE<dim>::computeRHS(const MatrixFree<dim,double> &data, 
				     std::vector<vectorType*> &dst, 
				     const std::vector<vectorType*> &src,
				     const std::pair<unsigned int,unsigned int> &cell_range) const{
}

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

//solve each time increment
template <int dim>
void MatrixFreePDE<dim>::solveIncrement(){
  //log time
  computing_timer.enter_section("matrixFreePDE: solveIncrements");
  Timer time; 

  //updateRHS
  updateRHS();
  
  //solve for each field
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    
    //Parabolic (first order derivatives in time) fields
    if (fields[fieldIndex].pdetype==PARABOLIC){
      //explicit-time step each DOF
      for (unsigned int dof=0; dof<solutionSet[fieldIndex]->local_size(); ++dof){
	solutionSet[fieldIndex]->local_element(dof)+=\
	  invM.local_element(dof)*timeStep*residualSet[fieldIndex]->local_element(dof);
      }
      char buffer[200];
      sprintf(buffer, "field '%s' [explicit solve]: current solution: %12.6e, current residual:%12.6e\n", \
	      fields[fieldIndex].name.c_str(),				\
	      solutionSet[fieldIndex]->l2_norm(),			\
	      residualSet[fieldIndex]->l2_norm()); 
      pcout<<buffer; 
    }
    
    //Elliptic (time-independent) fields
    else if (fields[fieldIndex].pdetype==ELLIPTIC){
      //implicit solve
      SolverControl solver_control(maxSolverIterations, relSolverTolerance*residualSet[fieldIndex]->l2_norm());
      solverType<vectorType> solver(solver_control);
      try{
	//solver.solve(*this, *solutionSet[fieldIndex], *residualSet[fieldIndex], IdentityMatrix(solutionSet[fieldIndex]->size()));
      }
      catch (...) {
	pcout << "\nWarning: solver did not converge as per set tolerances. consider increasing maxSolverIterations or decreasing relSolverTolerance.\n";
      }
      char buffer[200];
      sprintf(buffer, "field '%s' [implicit solve]: initial residual:%12.6e, current residual:%12.6e, nsteps:%u, tolerance criterion:%12.6e\n",\
	      fields[fieldIndex].name.c_str(),				\
	      solver_control.initial_value(),				\
	      solver_control.last_value(),				\
	      solver_control.last_step(), solver_control.tolerance()); 
      pcout<<buffer; 
    }
    
    //Hyperbolic (second order derivatives in time) fields and general
    //non-linear PDE types not yet implemented
    else{
      pcout << "matrixFreePDE.h: unknown field pdetype\n";
      exit(-1);
    }
  }
  pcout << "wall time: " << time.wall_time() << "s\n";
  //log time
  computing_timer.exit_section("matrixFreePDE: solveIncrements"); 
}

//output results
template <int dim>
void MatrixFreePDE<dim>::outputResults(){
  //log time
  computing_timer.enter_section("matrixFreePDE: output");
  
  //create DataOut object
  DataOut<dim> data_out;

  //loop over fields
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    //apply constraints
    constraintsSet[fieldIndex]->distribute (*solutionSet[fieldIndex]);
    //sync ghost DOF's
    solutionSet[fieldIndex]->update_ghost_values();
    //mark field as scalar/vector
    std::vector<DataComponentInterpretation::DataComponentInterpretation> dataType \
      (fields[fieldIndex].numComponents,				\
       (fields[fieldIndex].type==SCALAR ?				\
	DataComponentInterpretation::component_is_scalar:		\
	DataComponentInterpretation::component_is_part_of_vector));
    //add field to data_out
    data_out.add_data_vector(*dofHandlersSet[fieldIndex], *solutionSet[fieldIndex], fields[fieldIndex].name.c_str(), dataType);  
  }
  data_out.build_patches ();

  //write to results file
  //file name
  const std::string filename = "solution-" + \
    Utilities::int_to_string (currentIncrement, std::ceil(std::log10(totalIncrements))+1);
  //create file stream
  std::ofstream output ((filename +					\
			 "." + Utilities::int_to_string (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), \
							 std::ceil(std::log10(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)))+1) \
			 + ".vtu").c_str());
  //write to file
  pcout << filename.c_str() << std::endl;
  data_out.write_vtu (output);
  //create pvtu record
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0){
    std::vector<std::string> filenames;
    for (unsigned int i=0;i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i)
      filenames.push_back ("solution-" +				\
			   Utilities::int_to_string (currentIncrement, std::ceil(std::log10(totalIncrements))+1) \
			   + "." +					\
			   Utilities::int_to_string (i, std::ceil(std::log10(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)))+1) \
			   + ".vtu");
    std::ofstream master_output ((filename + ".pvtu").c_str());
    data_out.write_pvtu_record (master_output, filenames);
  }
  pcout << "Output written to: " << (filename + ".pvtu").c_str() << "\n\n";
  
  //log time
  computing_timer.exit_section("matrixFreePDE: output"); 
}

//solve BVP
template <int dim>
void MatrixFreePDE<dim>::solve(){
  //log time
  computing_timer.enter_section("matrixFreePDE: solve"); 
  
  //time dependent BVP
  if (isTimeDependentBVP){
    //initialize time step variables
    timeStep=timeStepV;
    finalTime=finalTimeV;
    totalIncrements=totalIncrementsV;    
    //output initial conditions for time dependent BVP
    if (writeOutput) outputResults();

    //time step
    for (currentIncrement=1; currentIncrement<totalIncrements; ++currentIncrement){
      //increment current time
      currentTime+=timeStep;
      pcout << "\ntime increment:" << currentIncrement << "  time: " << currentTime << "\n";
      if (currentTime>=finalTime){
	pcout << "\ncurrentTime>=finalTime. Ending time stepping\n";
	break;
      }
      //solve time increment
      solveIncrement();
      //output results to file
      if ((writeOutput) && (currentIncrement%skipOutputSteps==0)){
	outputResults();
      }
    }
  }
  //time independent BVP
  else{
    //solve
    solveIncrement();
    //output results to file
    if ((writeOutput) && (currentIncrement%skipOutputSteps==0)){
      outputResults();
    }
  }

  //log time
  computing_timer.exit_section("matrixFreePDE: solve"); 
}

#endif
