//computeInvM() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"
#include <numeric>

//compute inverse of the diagonal mass matrix and store in vector invM
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::computeInvM(){
	//initialize  invM
	bool invMInitialized=false;
	unsigned int parabolicFieldIndex=0;
	for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
		if (fields[fieldIndex].pdetype==EXPLICIT_TIME_DEPENDENT || fields[fieldIndex].pdetype==AUXILIARY){
			matrixFreeObject.initialize_dof_vector (invM, fieldIndex);
			parabolicFieldIndex=fieldIndex;
			invMInitialized=true;
			break;
		}
	}
	//check if invM initialized
	if (!invMInitialized){
		pcout << "matrixFreePDE.h: no PARABOLIC field... hence setting parabolicFieldIndex to 0 and marching ahead withn invM computation\n";
		//exit(-1);
	}

	//compute invM
	matrixFreeObject.initialize_dof_vector (invM, parabolicFieldIndex);
	invM=0.0;

	//select gauss lobatto quadrature points which are suboptimal but give diagonal M
	if (fields[parabolicFieldIndex].type==SCALAR){
		VectorizedArray<double> one = make_vectorized_array (1.0);
		FEEvaluation<dim,degree> fe_eval(matrixFreeObject, parabolicFieldIndex);
		const unsigned int n_q_points = fe_eval.n_q_points;
		for (unsigned int cell=0; cell<matrixFreeObject.n_cell_batches(); ++cell){
			fe_eval.reinit(cell);
			for (unsigned int q=0; q<n_q_points; ++q){
				fe_eval.submit_value(one,q);
			}
			fe_eval.integrate (true,false);
			fe_eval.distribute_local_to_global (invM);
		}
	}
	else {
		dealii::Tensor<1, dim, dealii::VectorizedArray<double> > oneV;
		for (unsigned int i=0;i<dim;i++){
			oneV[i] = 1.0;
		}

		FEEvaluation<dim,degree,degree+1,dim> fe_eval(matrixFreeObject, parabolicFieldIndex);

		const unsigned int n_q_points = fe_eval.n_q_points;
		for (unsigned int cell=0; cell<matrixFreeObject.n_cell_batches(); ++cell){
			fe_eval.reinit(cell);
			for (unsigned int q=0; q<n_q_points; ++q){
				fe_eval.submit_value(oneV,q);
			}
			fe_eval.integrate (true,false);
			fe_eval.distribute_local_to_global (invM);
		}
	}


	invM.compress(VectorOperation::add);

    // Calculate the volume of the smallest cell to prevent a non-zero value of invM being
    // confused for a near zero value (which can happen if the domain size is 1e-6 or below)
    std::vector<double> min_element_length;
    for (unsigned int d=0; d<dim; d++){
        int num_elements = userInputs.subdivisions.at(d)*dealii::Utilities::fixed_power<2>(userInputs.max_refinement_level);
        min_element_length.push_back(userInputs.domain_size[d]/double(num_elements));
    }
    double min_cell_volume = std::accumulate(begin(min_element_length), end(min_element_length), 1, std::multiplies<double>());

	//invert mass matrix diagonal elements
#if (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR < 4)
        for (unsigned int k=0; k<invM.local_size(); ++k){
#else
	for (unsigned int k=0; k<invM.locally_owned_size(); ++k){
#endif
		if (std::abs(invM.local_element(k))>1.0e-15 * min_cell_volume){
			invM.local_element(k) = 1./invM.local_element(k);
		}
		else{
			invM.local_element(k) = 0;
		}
	}
	pcout << "computed mass matrix (using FE space for field: " << parabolicFieldIndex << ")\n";
}

#include "../../include/matrixFreePDE_template_instantiations.h"
