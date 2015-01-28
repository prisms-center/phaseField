//Matrix Free implementation of Precipitate evolution (Coupled CH+AC+Mechanics) equations
#ifndef COUPLEDCH_AC_MECHANICS_H
#define COUPLEDCH_AC_MECHANICS_H

//general headers
#include <fstream>
#include <sstream>

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

//initial condition functions
//Concentration initial conditions
template <int dim>
class InitialConditionC : public Function<dim>
{
public:
  InitialConditionC () : Function<dim>(1) {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  virtual double value (const Point<dim> &p, const unsigned int component = 0) const;
};

//Structural order parameter initial conditions
template <int dim>
class InitialConditionN : public Function<dim>
{
public:
  unsigned int index;
  InitialConditionN (const unsigned int _index) : Function<dim>(1), index(_index) 
  {
    std::srand(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)+1);
  }
  virtual double value (const Point<dim> &p, const unsigned int component = 0) const;
};

//
//Precipitate evolution Class
//
template <int dim>
class PrecipitateProblem:  public Subscriptor
{
public:
  PrecipitateProblem ();
  void run ();
private:
  void setup_system ();
  void solve ();
  void cgSolve(vectorType &x, const vectorType &b);
  void output_results (const unsigned int cycle);
  parallel::distributed::Triangulation<dim>      triangulation;
  FE_Q<dim>                        fe;
  FESystem<dim>                    feU;
  DoFHandler<dim>                  dof_handler, dof_handlerU;
  ConstraintMatrix                 constraints, constraintsU;
  IndexSet                         locally_relevant_dofs, locally_relevant_dofsU;
  vectorType                       U, residualU, R, P, H;
  vectorType                       C, residualC;
  std::vector<vectorType>          N, residualN;
  std::vector<vectorType*>         solutions, residuals;
  double                           setup_time;
  unsigned int                     increment;
  ConditionalOStream               pcout;

  //elasticity matrix
  Table<2, double> CIJ;

  //matrix free objects
  MatrixFree<dim,double>      data;
  vectorType invM;
  void setup_matrixfree ();
  void vmult (vectorType &dst, const vectorType &src) const;
  void updateRHS();
  void computeRHS(const MatrixFree<dim,double> &data, 
		  std::vector<vectorType*> &dst, 
		  const std::vector<vectorType*> &src,
		  const std::pair<unsigned int,unsigned int> &cell_range) const;
  void computeAx(const MatrixFree<dim,double> &data, 
		 vectorType &dst, 
		 const vectorType &src,
		 const std::pair<unsigned int,unsigned int> &cell_range) const;
};


//constructor
template <int dim>
PrecipitateProblem<dim>::PrecipitateProblem ()
  :
  Subscriptor(),
  triangulation (MPI_COMM_WORLD),
  fe (QGaussLobatto<1>(finiteElementDegree+1)),
  feU (FE_Q<dim>(QGaussLobatto<1>(finiteElementDegree+1)),dim),
  dof_handler (triangulation),
  dof_handlerU (triangulation),
  N(numStructuralOrderParameters), 
  residualN(numStructuralOrderParameters),
  pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
  CIJ(2*dim-1+dim/3,2*dim-1+dim/3)
{
  //get elasticity matrix
  double materialConstants[]=MaterialConstantsv;
  getCIJMatrix<dim>(MaterialModelv, materialConstants, CIJ, pcout);
  //add solution vectors to output array
  solutions.push_back(&C); residuals.push_back(&residualC);
  for (unsigned int i=0; i<numStructuralOrderParameters; i++){
    solutions.push_back(&N[i]); residuals.push_back(&residualN[i]);
  }
  solutions.push_back(&U); residuals.push_back(&residualU);
}

//right hand side vector computation for explicit time stepping
template <int dim>
void PrecipitateProblem<dim>::computeRHS (const MatrixFree<dim,double>  &data,
					  std::vector<vectorType*> &dst, 
					  const std::vector<vectorType*> &src,
					  const std::pair<unsigned int,unsigned int> &cell_range) const
{
  double Mn[numStructuralOrderParameters]=MnVals;
  double Kn0[3][3]=Kn0Tensor;
  double Kn1[3][3]=Kn1Tensor;
  double Kn2[3][3]=Kn2Tensor;

  //initialize vals vectors for C, N
  std::vector<FEEvaluation<dim,finiteElementDegree,finiteElementDegree+1,1,double>*> vals;
  for (unsigned int i=0; i<numStructuralOrderParameters+1; i++){
    vals.push_back(new FEEvaluation<dim,finiteElementDegree,finiteElementDegree+1,1,double>(data,0));
  }
  //initialize vals vectors for U 
  FEEvaluation<dim,finiteElementDegree,finiteElementDegree+1,dim,double> valU(data,1);
  
  //initialize constants
  VectorizedArray<double> constCx, constN, constNx;
  constCx=-Mc; constN=-Mn[0]; constNx=-Mn[0];
  double sf0Strain[3][3]=sf0StrainV;
  double sf1Strain[3][3]=sf1StrainV;
  double sf2Strain[3][3]=sf2StrainV;

  //loop over all "cells"
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //read C, N values for this cell
    for (unsigned int i=0; i<numStructuralOrderParameters+1; i++){
      vals[i]->reinit (cell);
      vals[i]->read_dof_values(*src[i]);
      vals[i]->evaluate (true,true,false);
    }
    //read U values for this cell
    valU.reinit (cell); 
    valU.read_dof_values(*src[numStructuralOrderParameters+1]); 
    valU.evaluate (false,true,false); 
    
    //loop over quadrature points
    for (unsigned int q=0; q<vals[0]->n_q_points; ++q){
      //fill c,n values
      VectorizedArray<double> c  = vals[0]->get_value(q), n[numStructuralOrderParameters];
      Tensor<1, dim, VectorizedArray<double> >  cx = vals[0]->get_gradient(q), nx[numStructuralOrderParameters];
      for (unsigned int i=0; i<numStructuralOrderParameters; i++){
	n[i]=vals[i+1]->get_value(q);
	nx[i]=vals[i+1]->get_gradient(q);
      }
      //fill u values
      Tensor<1, dim, gradType> ux = valU.get_gradient(q);     
      //submit u values
      Tensor<1, dim, gradType> Cux;
      
      //compute stress
      dealii::VectorizedArray<double> E[dim][dim], R[dim][dim];
      //compute cumulative stress free strain
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  E[i][j]= sf0Strain[i][j]*h0V+sf1Strain[i][j]*h1V+sf2Strain[i][j]*h2V;
	}
      }
      //compute C*E0
      computeStress<dim>(CIJ, E, R);
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  Cux[1,i][1,j] = R[i][j];
	}
      }
      valU.submit_gradient(Cux,q);
      
      //C*(E-E0)*(sfStrain)
      //compute E2=E-E0
      dealii::VectorizedArray<double> E2[dim][dim];
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  E2[i][j]= 0.5*(ux[i][j]+ux[j][i])-(sf0Strain[i][j]*h0V+sf1Strain[i][j]*h1V+sf2Strain[i][j]*h2V);
	}
      }
      //compute R=C*(E-E0)
      computeStress<dim>(CIJ, E2, R);
      VectorizedArray<double> CEE0=make_vectorized_array(0.0);
      VectorizedArray<double> CEE1=make_vectorized_array(0.0);
      VectorizedArray<double> CEE2=make_vectorized_array(0.0);
      //compute R*sfStrain
      for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
	  CEE0+=R[i][j]*sf0Strain[i][j];
	  CEE1+=R[i][j]*sf1Strain[i][j];
	  CEE2+=R[i][j]*sf2Strain[i][j];	
	}
      }
        
      //submit c, n values
      for (unsigned int i=0; i<numStructuralOrderParameters+1; i++){
	if (i==0) vals[i]->submit_gradient(constCx*rcxV,q);
	else{
	  gradType rnx;
	  for (unsigned int a=0; a<dim; a++) {
	    rnx[1,a]=make_vectorized_array(0.0);
	    for (unsigned int b=0; b<dim; b++)
	      if (i==1) rnx[1,a]+=Kn0[a][b]*rn0xV[1,b];
	      else if (i==2) rnx[1,a]+=Kn1[a][b]*rn1xV[1,b];
	      else if (i==3) rnx[1,a]+=Kn2[a][b]*rn2xV[1,b];
	  }
	  //Assemble rnV-C*(E-E0)*(sfStrain)*hnV to value
	  if (i==1) vals[i]->submit_value(constN*(rn0V-CEE0*h0nV),q);
	  else if (i==2) vals[i]->submit_value(constN*(rn1V-CEE1*h1nV),q);
	  else if (i==3) vals[i]->submit_value(constN*(rn2V-CEE2*h2nV),q);
	  vals[i]->submit_gradient(constNx*rnx,q);
	}
      }  
    }
    //integrate c, n values
    for (unsigned int i=0; i<numStructuralOrderParameters+1; i++){
      if (i==0) vals[i]->integrate(false,true);
      else vals[i]->integrate(true,true);
      vals[i]->distribute_local_to_global (*dst[i]); 
    }
    //integrate u values
    valU.integrate(false, true); 
    valU.distribute_local_to_global(*dst[numStructuralOrderParameters+1]);
  }
  //release memory
  for (unsigned int i=0; i<numStructuralOrderParameters+1; i++){
    delete[]  vals[i];
  }
}

//Matrix-vector product for mechanics problem
template <int dim>
void PrecipitateProblem<dim>::computeAx (const MatrixFree<dim,double>  &data,
					 vectorType &dst, 
					 const vectorType &src,
					 const std::pair<unsigned int,unsigned int> &cell_range) const
{
  //initialize vals vectors
  FEEvaluation<dim,finiteElementDegree,finiteElementDegree+1,dim,double> vals(data,1);
  
  //loop over all "cells"
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    //read U values for this cell
    vals.reinit (cell); 
    vals.read_dof_values(src); 
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
	  Cux[1,i][1,j] = R[i][j];
	}
      }
      //compute residuals
      vals.submit_gradient(Cux,q);
    }
    vals.integrate(false, true); 
    vals.distribute_local_to_global(dst);
  }
}

// Matrix free data structure vmult operations.
template <int dim>
void PrecipitateProblem<dim>::vmult (vectorType &dst, const vectorType &src) const{
  dst=0.0;
  data.cell_loop (&PrecipitateProblem::computeAx, this, dst, src);
}

template <int dim>
void PrecipitateProblem<dim>::updateRHS (){
  for (unsigned int i=0; i<residuals.size(); i++)
    (*residuals[i])=0.0;
  data.cell_loop (&PrecipitateProblem::computeRHS, this, residuals, solutions);
}

//setup matrixfree data structures
template <int dim>
void PrecipitateProblem<dim>::setup_matrixfree (){
  //setup the matrix free object
  typename MatrixFree<dim,double>::AdditionalData additional_data;
  additional_data.mpi_communicator = MPI_COMM_WORLD;
  additional_data.tasks_parallel_scheme = MatrixFree<dim,double>::AdditionalData::partition_partition;
  additional_data.mapping_update_flags = (update_gradients | update_JxW_values);
  QGaussLobatto<1> quadrature (finiteElementDegree+1);
  //
  std::vector<const DoFHandler<dim>*> dof_handler_vector;
  dof_handler_vector.push_back(&dof_handler);
  dof_handler_vector.push_back(&dof_handlerU);
  //
  std::vector<const ConstraintMatrix*> constraint_vector;
  constraint_vector.push_back(&constraints);
  constraint_vector.push_back(&constraintsU);
  //
  data.reinit (dof_handler_vector, constraint_vector, quadrature, additional_data);

  //compute  invM
  data.initialize_dof_vector (invM, 0); 
  VectorizedArray<double> one = make_vectorized_array (1.0);
  
  //select gauss lobatto quad points which are suboptimal but give diogonal M 
  FEEvaluationGL<dim,finiteElementDegree,1,double> fe_eval(data,0);
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
void PrecipitateProblem<dim>::setup_system ()
{
  Timer time;
  time.start ();
  setup_time = 0;
  dof_handler.distribute_dofs (fe);   dof_handlerU.distribute_dofs (feU);
  pcout << "Number of global active cells: " << triangulation.n_global_active_cells() << std::endl;
  pcout << "Number of degrees of freedom:  " << (dof_handler.n_dofs() + dof_handlerU.n_dofs()) << std::endl;
  //setup dofhandler for C,N
  DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
  constraints.clear();
  constraints.reinit (locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  constraints.close();
  //setup dofhandler for U
  DoFTools::extract_locally_relevant_dofs (dof_handlerU, locally_relevant_dofsU);
  constraintsU.clear();
  constraintsU.reinit (locally_relevant_dofsU);
  DoFTools::make_hanging_node_constraints (dof_handlerU, constraintsU);
  VectorTools::interpolate_boundary_values (dof_handlerU, 0, ZeroFunction<dim>(dim), constraintsU);
  constraintsU.close();
  
  setup_time += time.wall_time();
  pcout << "Distribute DoFs & B.C.      "
	<<  time.wall_time() << "s" << std::endl;
  time.restart();

  //initialize matrix free data structure
  setup_matrixfree();

  //data structures
  data.initialize_dof_vector(C,0);
  residualC.reinit(C);
  for (unsigned int i=0; i<numStructuralOrderParameters; i++){
    N[i].reinit(C); residualN[i].reinit(C);
  }
  //data structures for U 
  data.initialize_dof_vector(U,1);
  residualU.reinit(U); 
  R.reinit(U); 
  P.reinit(U);
  H.reinit(U);  

  //initial Condition
  VectorTools::interpolate (dof_handler,InitialConditionC<dim> (),C);
  for (unsigned int i=0; i<numStructuralOrderParameters; i++){
    VectorTools::interpolate (dof_handler,InitialConditionN<dim> (i),N[i]);
  }
  U=0.0;

  //timing
  setup_time += time.wall_time();
  pcout << "Wall time for matrix-free setup:    "
	<< time.wall_time() << "s" << std::endl;
}

//CG solver for the mechanics problem
template <int dim>
void PrecipitateProblem<dim>::cgSolve(vectorType &x, const vectorType &b){
  char buffer[250];
  Timer time; double t=0, t1;
  double rtol=relativeResTol, rtolAbs=absoluteResTol, utol=dUTol;
  unsigned int maxIterations=4000, iterations=0;
  double res, res0, resOld;
  if (!x.all_zero()){
    //R=Ax-b
    t1=time.wall_time();
    vmult(R,x);
    t+=time.wall_time()-t1;
    R.add(-1.,b);
  }
  else
    R.equ(-1.,b); //R=-b
  P.equ(-1.,R); //P=-R=b-Ax
  res0 = res = R.l2_norm();
  //check convergence of res
  //sprintf(buffer,"InitialRes:%12.6e\n", res); pcout<<buffer;
  if (res<rtolAbs){sprintf(buffer,"Converged in absolute tolerance: initialRes:%12.6e, res:%12.6e, absTolCriterion:%12.6e, incs:%u\n", res, res, rtolAbs, 0); pcout<<buffer; return;}
  double alpha, beta;
  while(iterations<maxIterations){
    //compute alpha
    t1=time.wall_time();
    vmult(H,P);
    t+=time.wall_time()-t1;
    alpha=P*H;
    alpha=(res*res)/alpha;
    //compute R_(k+1), X_(k+1)
    R.add(alpha,H); //R=R+alpha*A*P
    x.add(alpha,P); //X=X+alpha*P
    resOld=res;
    res = R.l2_norm();
    //sprintf(buffer, "%12.6e ", alpha*P.linfty_norm()); pcout<<buffer;
    if (res<=rtolAbs){sprintf(buffer, "Elasticity equations converged in absolute R tolerance: initialRes:%12.6e, currentRes:%12.6e, tol criterion:%12.6e, incs:%u\n", res, res, rtolAbs, iterations); pcout<<buffer; sprintf(buffer, "time in mat-mult:%12.6e\n", t); pcout<<buffer; return;}
    else if ((res/res0)<=rtol){printf("Elasticity equations converged in relative R tolerance: initialRes:%12.6e, currentRes:%12.6e, relTol:%12.6e, tol criterion:%12.6e, dU:%12.6e, incs:%u\n", res0, res, res/res0, rtol, alpha*P.linfty_norm(), iterations); printf("time in mat-mult:%12.6e\n", t); return;}
    else if ((alpha*P.linfty_norm())<utol){sprintf(buffer, "Elasticity equations converged in dU tolerance: initialRes:%12.6e, currentRes:%12.6e, relTol:%12.6e, tol criterion:%12.6e, dU:%12.6e, incs:%u\n", res0, res, res/res0, utol, alpha*P.linfty_norm(), iterations); pcout<<buffer; sprintf(buffer, "time in mat-mult: %12.6es\n", t); pcout<<buffer; return;}
    //compute beta
    beta = (res*res)/(resOld*resOld);
    P.sadd(beta,-1.0,R);
    iterations++;
  }
  sprintf(buffer,"Res:%12.6e, incs:%u\n", res, iterations); pcout<<buffer;
}


//solve
template <int dim>
void PrecipitateProblem<dim>::solve ()
{
  Timer time; 
  updateRHS();
  //compute c=c0-mobility*dt*f'(c0)
  for (unsigned int j=0; j<C.local_size(); ++j){
    C.local_element(j)+=invM.local_element(j)*dt*residualC.local_element(j);
  }
  //compute n0-mobility*dt*f'(n0)
  for (unsigned int i=0; i<numStructuralOrderParameters; i++){
    for (unsigned int j=0; j<N[i].local_size(); ++j){
      N[i].local_element(j)+=invM.local_element(j)*dt*residualN[i].local_element(j);
    }
  }
  //compute u
  if (increment%skipElasticitySteps==0) {
    cgSolve(U,residualU); 
  }
  pcout << "solve wall time: " << time.wall_time() << "s\n";
}
  
//output 
template <int dim>
void PrecipitateProblem<dim>::output_results (const unsigned int cycle)
{
  //apply constraints
  constraints.distribute (C); C.update_ghost_values();
  for (unsigned int i=0; i<numStructuralOrderParameters; i++){
    constraints.distribute (N[i]);
    N[i].update_ghost_values();
  }
  constraintsU.distribute (U);
  U.update_ghost_values();
  //Create DataOut object
  DataOut<dim> data_out;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> dataType(1, DataComponentInterpretation::component_is_scalar);
  //data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (dof_handler, C, "c", dataType);
  for (unsigned int i=0; i<numStructuralOrderParameters; i++){
    const std::string fieldname = "eta" + Utilities::int_to_string (i+1, 1);
    data_out.add_data_vector (dof_handler, N[i], fieldname.c_str(), dataType);
  }
  std::vector<std::string> solution_names (dim, "u");
  std::vector<DataComponentInterpretation::DataComponentInterpretation> dataTypeU(dim, DataComponentInterpretation::component_is_part_of_vector);
  data_out.add_data_vector (dof_handlerU, U, solution_names, dataTypeU); 
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
  pcout << "Output written to:" << filename.c_str() << "\n\n";
}


//run
template <int dim>
void PrecipitateProblem<dim>::run ()
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
      pcout << "\nTime increment:" << increment << " T: " << dt*increment << ", Wall time: " << time.wall_time() << "s\n";
      //call solve method
      solve ();
      //write results to file
      if ((writeOutput) && (increment%skipOutputSteps==0)) output_results (increment);
    };
  pcout << "Total  time    "  << time.wall_time() << "s" << std::endl; 
  pcout << "PerInc time    "  << time.wall_time()/numIncrements << "s" << std::endl; 
}

#endif
