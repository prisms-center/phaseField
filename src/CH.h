//Matrix Free implementation of Cahn-Hilliard spinodal decomposition mixed (split) formulation 
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

//
//Cahn-Hilliard class
//
template <int dim>
class CahnHilliardProblem:  public Subscriptor
{
public:
  CahnHilliardProblem ();
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
  vectorType                       Mu, residualMu;
  vectorType                       C,  residualC;
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
  bool updateMuRHS;
};


//constructor
template <int dim>
CahnHilliardProblem<dim>::CahnHilliardProblem ()
  :
  Subscriptor(),
  triangulation (MPI_COMM_WORLD),
  fe (QGaussLobatto<1>(finiteElementDegree+1)),
  dof_handler (triangulation),
  pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
{
  solutions.push_back(&Mu); residuals.push_back(&residualMu);
  solutions.push_back(&C); residuals.push_back(&residualC);
}

//right hand side vector computation for explicit time stepping
template <int dim>
void CahnHilliardProblem<dim>::computeRHS (const MatrixFree<dim,double>  &data,
					  std::vector<vectorType*> &dst, 
					  const std::vector<vectorType*> &src,
					  const std::pair<unsigned int,unsigned int> &cell_range) const
{
  //initialize vals vectors
  FEEvaluation<dim,finiteElementDegree> valMu(data), valC(data);
  VectorizedArray<double> constMux, constCx;
  constMux=Kc; constCx=-Mc*dt;
  
  //loop over all "cells"
  for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){
    
    //read Mu, C values for this cell
    valMu.reinit (cell); valC.reinit (cell);
    valMu.read_dof_values(*src[0]); valMu.evaluate (true,true,false); 
    valC.read_dof_values(*src[1]);  valC.evaluate (true,true,false); 
    
    //loop over quadrature points
    for (unsigned int q=0; q<valMu.n_q_points; ++q){
      VectorizedArray<double> mu  = valMu.get_value(q), c  = valC.get_value(q);
      Tensor<1, dim, VectorizedArray<double> >  mux = valMu.get_gradient(q), cx = valC.get_gradient(q);
      
      //compute residuals
      if (updateMuRHS){
	valMu.submit_value(rmuV,q);
	valMu.submit_gradient(constMux*rmuxV,q);
      }
      else{
	valC.submit_value(rcV,q);
	valC.submit_gradient(constCx*rcxV,q);
      }
    }
    if (updateMuRHS){
      valMu.integrate(true, true); valMu.distribute_local_to_global(*dst[0]);
    }
    else{
      valC.integrate(true, true);  valC.distribute_local_to_global(*dst[1]);
    }
  }
}

template <int dim>
void CahnHilliardProblem<dim>::updateRHS (){
  //initialize residuals to zero
  for (unsigned int i=0; i<residuals.size(); i++)
    (*residuals[i])=0.0;
  //loop over all cells to compute residuals
  data.cell_loop (&CahnHilliardProblem::computeRHS, this, residuals, solutions);
}

//setup matrixfree data structures
template <int dim>
void CahnHilliardProblem<dim>::setup_matrixfree (){
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
void CahnHilliardProblem<dim>::setup_system ()
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
  data.initialize_dof_vector(Mu);
  C.reinit(Mu);
  residualMu.reinit(Mu); 
  residualC.reinit(Mu);

  //initial Condition
  Mu=0.0;
  VectorTools::interpolate (dof_handler,InitialConditionC<dim> (), C);
  
  //timing
  setup_time += time.wall_time();
  pcout << "Wall time for matrix-free setup:    "
	<< time.wall_time() << "s" << std::endl;
}

//solve
template <int dim>
void CahnHilliardProblem<dim>::solve ()
{
  Timer time; 
  //compute C^(n+1)
  updateMuRHS=false; updateRHS(); 
  for (unsigned int j=0; j<C.local_size(); ++j)
    C.local_element(j)=invM.local_element(j)*residualC.local_element(j);
 
  //compute Mu^(n+1)
  updateMuRHS=true; updateRHS(); 
  for (unsigned int j=0; j<Mu.local_size(); ++j)
    Mu.local_element(j)=invM.local_element(j)*residualMu.local_element(j);

  pcout << "solve wall time: " << time.wall_time() << "s\n";
}
  
//output 
template <int dim>
void CahnHilliardProblem<dim>::output_results (const unsigned int cycle)
{
  constraints.distribute (Mu); 
  Mu.update_ghost_values();
  constraints.distribute (C);
  C.update_ghost_values();
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (Mu, "mu");  data_out.add_data_vector (C, "c");
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
void CahnHilliardProblem<dim>::run ()
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

  //compute Mu^0
  updateMuRHS=true; updateRHS(); 
  for (unsigned int j=0; j<Mu.local_size(); ++j)
    Mu.local_element(j)=invM.local_element(j)*residualMu.local_element(j);

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

