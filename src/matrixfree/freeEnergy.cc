//computeFreeEnergyValue() virtual method for MatrixFreePDE class

#ifndef COMPUTEFREEENERGY_MATRIXFREE_H
#define COMPUTEFREEENERGY_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//compute integrated free energy value over the domain
template <int dim>
void MatrixFreePDE<dim>::computeFreeEnergyValue(std::vector<double>& freeEnergyValues){
}

// output the integrated free energies into a text file
template <int dim>
void MatrixFreePDE<dim>::outputFreeEnergy(std::vector<double>& freeEnergyValues){

	  std::ofstream output_file("./freeEnergy.txt");
	  output_file.precision(10);
	  std::ostream_iterator<double> output_iterator(output_file, "\n");
	  std::copy(freeEnergyValues.begin(), freeEnergyValues.end(), output_iterator);
}


#endif





