//outputResults() method for MatrixFreePDE class

#ifndef OUTPUTRESULTS_MATRIXFREE_H
#define OUTPUTRESULTS_MATRIXFREE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//output results
template <int dim>
void MatrixFreePDE<dim>::outputResults(){
  //log time
  computing_timer.enter_section("matrixFreePDE: output");
  
  //create DataOut object
  DataOut<dim> data_out;

  //loop over fields
  for(unsigned int fieldIndex=0; fieldIndex<fields.size(); fieldIndex++){
    //apply constraints
    constraintsSet[fieldIndex]->distribute (*solutionSet[fieldIndex]);
    //sync ghost DOF's
    solutionSet[fieldIndex]->update_ghost_values();
    //mark field as scalar/vector
    std::vector<DataComponentInterpretation::DataComponentInterpretation> dataType \
      (fields[fieldIndex].numComponents,				\
       (fields[fieldIndex].type==SCALAR ?				\
	DataComponentInterpretation::component_is_scalar:		\
	DataComponentInterpretation::component_is_part_of_vector));
    //add field to data_out
    std::vector<std::string> solutionNames (fields[fieldIndex].numComponents, fields[fieldIndex].name.c_str());
    data_out.add_data_vector(*dofHandlersSet[fieldIndex], *solutionSet[fieldIndex], solutionNames, dataType);  
  }
  data_out.build_patches ();

  //write to results file
  //file name
  const std::string filename = "solution-" + \
    Utilities::int_to_string (currentIncrement, std::ceil(std::log10(totalIncrements))+1);
  //create file stream
  std::ofstream output ((filename +					\
			 "." + Utilities::int_to_string (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), \
							 std::ceil(std::log10(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)))+1) \
			 + ".vtu").c_str());
  //write to file
  data_out.write_vtu (output);
  //create pvtu record
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0){
    std::vector<std::string> filenames;
    for (unsigned int i=0;i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i)
      filenames.push_back ("solution-" +				\
			   Utilities::int_to_string (currentIncrement, std::ceil(std::log10(totalIncrements))+1) \
			   + "." +					\
			   Utilities::int_to_string (i, std::ceil(std::log10(Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD)))+1) \
			   + ".vtu");
    std::ofstream master_output ((filename + ".pvtu").c_str());
    data_out.write_pvtu_record (master_output, filenames);
  }
  pcout << "Output written to: " << (filename + ".pvtu").c_str() << "\n\n";
  
  //log time
  computing_timer.exit_section("matrixFreePDE: output"); 
}

#endif

