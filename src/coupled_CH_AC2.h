//Matrix Free implementation of Precipitate evolution (Coupled CH+AC) equations
//general headers
#include <fstream>
#include <sstream>

//Code specific initializations.  
typedef dealii::parallel::distributed::Vector<double> vectorType;
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
  InitialConditionN (const unsigned int component) : Function<dim>(1) 
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
  void output_results (const unsigned int cycle);
  parallel::distributed::Triangulation<dim>      triangulation;
  FE_Q<dim>                        fe;
  DoFHandler<dim>                  dof_handler;
  ConstraintMatrix                 constraints;
  IndexSet                         locally_relevant_dofs;
  vectorType                       C, residualC;
  vectorType                       N, residualN;
  std::vector<vectorType*>         solutions, residuals;
  double                           setup_time;
  unsigned int                     increment;
  ConditionalOStream               pcout;

  //matrix free objects
  MatrixFree<dim,double>      data;
  vectorType invM;
  void setup_matrixfree ();
  void updateRHS();
  void computeRHS(const MatrixFree<dim,double> &data, 
		  std::vector<vectorType*> &dst, 
		  const std::vector<vectorType*> &src,
		  const std::pair<unsigned int,unsigned int> &cell_range) const;
};


//constructor
template <int dim>
PrecipitateProblem<dim>::PrecipitateProblem ()
  :
  Subscriptor(),
  triangulation (MPI_COMM_WORLD),
  fe (QGaussLobatto<1>(finiteElementDegree+1)),
  dof_handler (triangulation),
  pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
{
  solutions.push_back(&C); residuals.push_back(&residualC);
  solutions.push_back(&N); residuals.push_back(&residualN);
}

//right hand side vector computation for explicit time stepping
template <int dim>
void PrecipitateProblem<dim>::computeRHS (const MatrixFree<dim,double>  &data,
					  std::vector<vectorType*> &dst, 
					  const std::vector<vectorType*> &src,
					  const std::pair<unsigned int,unsigned int> &cell_range) const
{
  //initialize vals vectors
  FEEvaluation<dim,finiteElementDegree> valC(data), valN(data);
  VectorizedArray<double> constCx, constN, constNx;
  constCx=-Mc; constN=-Mn; constNx=-Mn*Kn;
  //loop over all "cells"
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    valC.reinit (cell); valN.reinit (cell);
    valC.read_dof_values(*src[0]); valN.read_dof_values(*src[1]);
    valC.evaluate (true,true,false);  valN.evaluate (true,true,false);
    for (unsigned int q=0; q<valC.n_q_points; ++q){
      VectorizedArray<double> c  = valC.get_value(q), n=valN.get_value(q);
      Tensor<1, dim, VectorizedArray<double> >  cx = valC.get_gradient(q), nx = valN.get_gradient(q);
      valC.submit_gradient(constCx*rcxV,q);
      valN.submit_value(constN*rnV,q);
      valN.submit_gradient(constNx*rnxV,q);
    }
    valC.integrate(false,true);
    valN.integrate(true,true);
    valC.distribute_local_to_global (*dst[0]); 
    valN.distribute_local_to_global (*dst[1]); 
  }
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
  data.reinit (dof_handler, constraints, quadrature, additional_data);

  //compute  invM
  data.initialize_dof_vector (invM); 
  VectorizedArray<double> one = make_vectorized_array (1.0);
  
  //select gauss lobatto quad points which are suboptimal but give diogonal M 
  FEEvaluationGL<dim,finiteElementDegree,1,double> fe_eval(data);
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
  //system_matrix.clear();
  dof_handler.distribute_dofs (fe);
  pcout << "Number of global active cells: " << triangulation.n_global_active_cells() << std::endl;
  pcout << "Number of degrees of freedom:  " << dof_handler.n_dofs() << std::endl;
  DoFTools::extract_locally_relevant_dofs (dof_handler, locally_relevant_dofs);
  constraints.clear();
  constraints.reinit (locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (dof_handler, constraints);
  
  //FIXME: Check Dirichlet BC
  //VectorTools::interpolate_boundary_values (dof_handler, 0, ZeroFunction<dim>(), constraints);
  constraints.close();
  setup_time += time.wall_time();
  pcout << "Distribute DoFs & B.C.      "
	<<  time.wall_time() << "s" << std::endl;
  time.restart();

  //initialize matrix free data structure
  setup_matrixfree();

  //data structures
  data.initialize_dof_vector(C);
  residualC.reinit(C);
  N.reinit(C); 
  residualN.reinit(C);

  //initial Condition
  VectorTools::interpolate (dof_handler,InitialConditionC<dim> (),C);
  VectorTools::interpolate (dof_handler,InitialConditionN<dim> (0),N);
    
  //timing
  setup_time += time.wall_time();
  pcout << "Wall time for matrix-free setup:    "
	<< time.wall_time() << "s" << std::endl;
}

//solve
template <int dim>
void PrecipitateProblem<dim>::solve ()
{
  Timer time; 
  //compute n0-mobility*dt*f'(n0)
  updateRHS();
  for (unsigned int j=0; j<C.local_size(); ++j)
    C.local_element(j)+=invM.local_element(j)*dt*residualC.local_element(j);
  for (unsigned int j=0; j<N.local_size(); ++j)
    N.local_element(j)+=invM.local_element(j)*dt*residualN.local_element(j);
  pcout << "solve wall time: " << time.wall_time() << "s\n";
}
  
//output 
template <int dim>
void PrecipitateProblem<dim>::output_results (const unsigned int cycle)
{
  constraints.distribute (C);
  constraints.distribute (N);
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (C, "c");
  data_out.add_data_vector (N, "eta"); 
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

