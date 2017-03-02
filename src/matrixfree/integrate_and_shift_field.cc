// This method needs to be generalized or removed. It probably belongs elsewhere. Exists historically for singlePrecipitateKKS.

#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
void  MatrixFreePDE<dim,degree>::computeIntegral(double& integratedField) {
	QGauss<dim>  quadrature_formula(degree+1);
	FE_Q<dim> FE (QGaussLobatto<1>(degree+1));
	FEValues<dim> fe_values (FE, quadrature_formula, update_values | update_JxW_values | update_quadrature_points);
	const unsigned int   n_q_points    = quadrature_formula.size();
	std::vector<double> cVal(n_q_points);

	typename DoFHandler<dim>::active_cell_iterator cell= dofHandlersSet[0]->begin_active(), endc = dofHandlersSet[0]->end();

	double value = 0.0;

	unsigned int fieldIndex;
	fieldIndex=getFieldIndex("c");

	for (; cell!=endc; ++cell) {
		if (cell->is_locally_owned()){
			fe_values.reinit (cell);

			fe_values.get_function_values(*solutionSet[fieldIndex], cVal);

			for (unsigned int q=0; q<n_q_points; ++q){
	    		value+=(cVal[q])*fe_values.JxW(q);
	    	}
		  }
	  }

	  value=Utilities::MPI::sum(value, MPI_COMM_WORLD);

	  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0){
	  std::cout<<"Integrated field: "<<value<<std::endl;
	  }

	  integratedField = value;
}

template <int dim, int degree>
void  MatrixFreePDE<dim,degree>::shiftConcentration(){
  //default trivial implementation.
}

#include "../../include/matrixFreePDE_template_instantiations.h"
