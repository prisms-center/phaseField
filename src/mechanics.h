//Matrix Free implementation of infinitesimal strain mechanics
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
private:
  void setup_system ();
  void solve ();
  void cgSolve(vectorType &x, const vectorType &b);
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

  //matrix free objects
  MatrixFree<dim,double>      data;
  vectorType invM;
  void setup_matrixfree ();
  void vmult (vectorType &dst, const vectorType &src) const;
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
  pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
{}

//right hand side vector computation for explicit time stepping
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
    vals.read_dof_values(src); 
    vals.evaluate (false,true,false); 
    
    //loop over quadrature points
    for (unsigned int q=0; q<vals.n_q_points; ++q){
      Tensor<1, dim, gradType> ux = vals.get_gradient(q);
      Tensor<1, dim, gradType> Cux;
      //compute C_(ijkl)*u_(k,l). Using minor symmetry of C tensor, C_(ijkl)=C_(ijlk).
      for (unsigned int i=0; i<dim; i++)
	for (unsigned int j=0; j<dim; j++)
	  for (unsigned int k=0; k<dim; k++)
	    for (unsigned int l=0; l<dim; l++){
	      if (updateRHSValue) {	
		double value= 2.0*(0.5 - (double)(std::rand() % 100 )/100.0);
		Cux[1,i][1,j] += CijklV*((double) (k==l))*make_vectorized_array(-1.0e-3*value);
	      }
	      else{
		Cux[1,i][1,j] += CijklV*ux[1,k][1,l];
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
  
  //FIXME: Check Dirichlet BC
  VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(dim), constraints);
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

  pcout << "size:" << U.size() << "/n";
  //initial Condition
  U=0.0;
  
  //timing
  setup_time += time.wall_time();
  pcout << "Wall time for matrix-free setup:    "
	<< time.wall_time() << "s" << std::endl;
}

template <int dim>
void MechanicsProblem<dim>::cgSolve(vectorType &x, const vectorType &b){
  char buffer[250];
  Timer time; double t=0, t1;
  double rtol=1.0e-10, rtolAbs=1.0e-15, utol=1.0e-6;
  unsigned int maxIterations=1000, iterations=0;
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
    if (res<=rtolAbs){sprintf(buffer, "Converged in absolute R tolerance: initialRes:%12.6e, currentRes:%12.6e, absTolCriterion:%12.6e, incs:%u\n", res, res, rtolAbs, iterations); pcout<<buffer; sprintf(buffer, "time in mat-mult:%12.6e\n", t); pcout<<buffer; return;}
    else if ((res/res0)<=rtol){printf("Converged in relative R tolerance: initialRes:%12.6e, currentRes:%12.6e, relTol:%12.6e, relTolCriterion:%12.6e, dU:%12.6e, incs:%u\n", res0, res, res/res0, rtol, alpha*P.linfty_norm(), iterations); printf("time in mat-mult:%12.6e\n", t); return;}
    else if ((alpha*P.linfty_norm())<rtol){sprintf(buffer, "Converged in dU tolerance: initialRes:%12.6e, currentRes:%12.6e, relTol:%12.6e, relTolCriterion:%12.6e, dU:%12.6e, incs:%u\n", res0, res, res/res0, rtol, alpha*P.linfty_norm(), iterations); pcout<<buffer; sprintf(buffer, "time in mat-mult: %12.6es\n", t); pcout<<buffer; return;}
    //compute beta
    beta = (res*res)/(resOld*resOld);
    P.sadd(beta,-1.0,R);
    iterations++;
  }
  sprintf(buffer,"Res:%12.6e, incs:%u\n", res, iterations); pcout<<buffer;
}

//solve
template <int dim>
void MechanicsProblem<dim>::solve ()
{
  Timer time; 
  updateRHS(); 
  cgSolve(U,residualU); 
  pcout << "solve wall time: " << time.wall_time() << "s\n";
}
  
//output 
template <int dim>
void MechanicsProblem<dim>::output_results (const unsigned int cycle)
{
  constraints.distribute (U);
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (U, "u");
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
  pcout << "Output written to:" << filename.c_str() << "\n\n";
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
      pcout << "\nTime increment:" << increment  << ", Wall time: " << time.wall_time() << "s\n";
      //call solve method
      solve ();
      //write results to file
      if ((writeOutput) && (increment%skipOutputSteps==0)) output_results (increment);
    };
  pcout << "Total  time    "  << time.wall_time() << "s" << std::endl; 
  pcout << "PerInc time    "  << time.wall_time()/numIncrements << "s" << std::endl; 
}

