//#include "../../include/postprocessor.h"
#include "../../include/matrixFreePDE.h"
#include "../../include/list_of_CIJ.h"

template <int dim,int degree>
void MatrixFreePDE<dim,degree>::computePostProcessedFields(std::vector<vectorType*> &postProcessedSet) const {


	// Zero out the postProcessedSet
	for(unsigned int fieldIndex=0; fieldIndex<userInputs.pp_number_of_variables; fieldIndex++){
		vectorType *U;
		U=new vectorType;
		postProcessedSet.push_back(U);
		matrixFreeObject.initialize_dof_vector(*U,  0); *U=0;
	}

	//call to integrate and assemble
	matrixFreeObject.cell_loop (&MatrixFreePDE<dim,degree>::getPostProcessedFields, this, postProcessedSet, solutionSet);

}

template <int dim,int degree>
void MatrixFreePDE<dim,degree>::getPostProcessedFields(const dealii::MatrixFree<dim,double> &data,
		std::vector<vectorType*> &dst,
		const std::vector<vectorType*> &src,
		const std::pair<unsigned int,unsigned int> &cell_range) const {

	//initialize FEEvaulation objects
	variableContainer<dim,degree,dealii::VectorizedArray<double> > variable_list(data,userInputs.varInfoListRHS);
	variableContainer<dim,degree,dealii::VectorizedArray<double> > pp_variable_list(data,userInputs.pp_varInfoList);

	//loop over cells
	for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){

		// Initialize, read DOFs, and set evaulation flags for each variable
		variable_list.reinit_and_eval(src, cell);
		pp_variable_list.reinit(cell);

		unsigned int num_q_points = variable_list.get_num_q_points();

		//loop over quadrature points
		for (unsigned int q=0; q<num_q_points; ++q){
			variable_list.q_point = q;
			pp_variable_list.q_point = q;
			
			dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc = variable_list.get_q_point_location();

			// Calculate the residuals
			postProcessedFields(variable_list,pp_variable_list,q_point_loc);

		}

		pp_variable_list.integrate_and_distribute(dst);
	}
}

#include "../../include/matrixFreePDE_template_instantiations.h"
