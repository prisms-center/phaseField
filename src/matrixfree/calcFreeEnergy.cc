//computeEnergy() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

//update RHS of each field
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::computeEnergy(){
  //log time
  computing_timer.enter_section("matrixFreePDE: computeEnergy");

  //call to integrate and assemble
  energy=0.0;
  energy_components.clear();
  energy_components.push_back(0.0);
  energy_components.push_back(0.0);
  energy_components.push_back(0.0);

  matrixFreeObject.cell_loop (&MatrixFreePDE<dim,degree>::getEnergy, this, residualSet, solutionSet);

  //add across all processors
  energy=Utilities::MPI::sum(energy, MPI_COMM_WORLD);
  for (unsigned int i=0; i<3; i++){
	  energy_components[i]=Utilities::MPI::sum(energy_components[i], MPI_COMM_WORLD);
  }
  pcout << "Energy: " << energy << std::endl;
  pcout << "Energy Components: " << energy_components[0] << " " << energy_components[1] << " " << energy_components[2] << " " << std::endl;
  freeEnergyValues.push_back(energy);
  //end log
  computing_timer.exit_section("matrixFreePDE: computeEnergy");
}

template <int dim, int degree>
void  MatrixFreePDE<dim,degree>::getEnergy(const MatrixFree<dim,double> &data,
				    std::vector<vectorType*> &dst,
				    const std::vector<vectorType*> &src,
				    const std::pair<unsigned int,unsigned int> &cell_range) {

    variableContainer<dim,degree,dealii::VectorizedArray<double> > variable_list(data,userInputs.varInfoListRHS);

	//loop over cells
    for (unsigned int cell=cell_range.first; cell<cell_range.second; ++cell){

        // Initialize, read DOFs, and set evaulation flags for each variable
        variable_list.reinit_and_eval(src, cell);

        unsigned int num_q_points = variable_list.get_num_q_points();

        dealii::AlignedVector<dealii::VectorizedArray<double> > JxW(num_q_points);

        variable_list.get_JxW(JxW);

        //loop over quadrature points
        for (unsigned int q=0; q<num_q_points; ++q){
            variable_list.q_point = q;
            dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc = variable_list.get_q_point_location();

            // Calculate the energy density
            energyDensity(variable_list,JxW[q],q_point_loc);
        }
    }

}

// output the integrated free energies into a text file
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::outputFreeEnergy(const std::vector<double>& freeEnergyValues) const{

	  std::ofstream output_file("./freeEnergy.txt");
	  output_file.precision(10);
	  std::ostream_iterator<double> output_iterator(output_file, "\n");
	  std::copy(freeEnergyValues.begin(), freeEnergyValues.end(), output_iterator);
}

#include "../../include/matrixFreePDE_template_instantiations.h"
