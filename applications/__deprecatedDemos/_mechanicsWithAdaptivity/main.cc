//mechanics problem implementation
//general headers
#include <fstream>
#include <sstream>

//dealIIheaders
#include "../../../include/dealIIheaders.h"

//mechanics problem headers
#include "parameters.h"
#include "../../../src/models/mechanics/mechanics.h"
#include "ICs_and_BCs.h"


//adaptive refinement control
template <int dim>
void MechanicsProblem<dim>::adaptiveRefine(unsigned int currentIncrement){
  this->refineMesh(currentIncrement);
}

//adaptive refinement criterion
template <int dim>
void MechanicsProblem<dim>::adaptiveRefineCriterion(){
  //Custom defined estimation criterion
  QGauss<dim>  quadrature(finiteElementDegree+1);
  FEValues<dim> fe_values (*this->FESet[refinementDOF], quadrature, update_values);
  const unsigned int   num_quad_points = quadrature.size();
  typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandlersSet2[refinementDOF]->begin_active(), endc = this->dofHandlersSet2[refinementDOF]->end();
  unsigned int vertices_per_cell=GeometryInfo<dim>::vertices_per_cell;
  for (;cell!=endc; ++cell){
    if (cell->is_locally_owned()){
      bool mark_refine = false;
      for (unsigned int i=0; i<vertices_per_cell; ++i){
	unsigned int nodeID=cell->vertex_dof_index(i,0);
	Point<3> feNodeGlobalCoord = cell->vertex(i);
	Point <3> center(50,50,50);
	if (center.distance(feNodeGlobalCoord)<20.0) mark_refine=true;
	//if ((feNodeGlobalCoord[0]<10.0) || (feNodeGlobalCoord[1]<10.0)) mark_refine=true;
      }
      if (mark_refine == true){
	cell->set_refine_flag();
      }
      else {
	//cell->set_coarsen_flag();
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
      MechanicsProblem<problemDIM> problem;
      problem.fields.push_back(Field<problemDIM>(VECTOR,  ELLIPTIC, "u"));
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
