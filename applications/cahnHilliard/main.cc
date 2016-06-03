//Cahn-Hilliard spinodal decomposition implementation
//general headers
#include "../../include/dealIIheaders.h"

//Cahn-Hilliard problem headers
#include "parameters.h"
#include "residuals.h"
#include "../../src/models/diffusion/CH.h"
#include "ICs_and_BCs.h"

//adaptive refinement control
template <int dim>
void CahnHilliardProblem<dim>::adaptiveRefine(unsigned int currentIncrement){
  if ((currentIncrement>0) && (currentIncrement%1000==0)){
    this->refineMesh(currentIncrement);
  }
}

//adaptive refinement criterion
template <int dim>
void CahnHilliardProblem<dim>::adaptiveRefineCriterion(){
  //Custom defined estimation criterion
  QGauss<dim>  quadrature(finiteElementDegree+1);
  FEValues<dim> fe_values (*this->FESet[refinementDOF], quadrature, update_values);
  const unsigned int   num_quad_points = quadrature.size();
  typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandlersSet2[refinementDOF]->begin_active(), endc = this->dofHandlersSet2[refinementDOF]->end();
  for (;cell!=endc; ++cell){
    if (cell->is_locally_owned()){
      fe_values.reinit (cell);
      std::vector<double> errorOut(num_quad_points);
      fe_values.get_function_values(*this->solutionSet[refinementDOF], errorOut);
      bool mark_refine = false;
      for (unsigned int q_point=0; q_point<num_quad_points; ++q_point){
    	  if ((errorOut[q_point]>0.11) && (errorOut[q_point]<0.99)){
    		  mark_refine = true;
    		  break;
    	  }
      }
      if (mark_refine == true){
    	  cell->set_refine_flag();
      }
      else {
    	  cell->set_coarsen_flag();
      }
    }
  }
}


//main
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,numbers::invalid_unsigned_int);
  try
    {
      deallog.depth_console(0);
      CahnHilliardProblem<problemDIM> problem;
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "mu"));
      problem.fields.push_back(Field<problemDIM>(SCALAR, PARABOLIC, "c"));
      problem.init (); 
      problem.solve();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  
  return 0;
}
