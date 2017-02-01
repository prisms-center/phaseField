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
  
  data_out.build_patches (finiteElementDegree);
  
  //write to results file
  //file name
  std::ostringstream cycleAsString;
  cycleAsString << std::setw(std::ceil(std::log10(totalIncrements))+1) << std::setfill('0') << currentIncrement;
  char vtuFileName[100], pvtuFileName[100];
  sprintf(vtuFileName, "solution-%s.%u.%s", cycleAsString.str().c_str(),Utilities::MPI::this_mpi_process(MPI_COMM_WORLD),outputFileType);
  sprintf(pvtuFileName, "solution-%s.p%s", cycleAsString.str().c_str(),outputFileType);
  std::ofstream output (vtuFileName);

  //write to file
  if (outputFileType == "vtu"){
	  data_out.write_vtu (output);
  }
  else if (outputFileType == "vtk"){
	  data_out.write_vtk (output);
  }
  else {
	  std::cout << "PRISMS-PF Error: The parameter 'outputFileType' must be either \"vtu\" or \"vtk\"" << std::endl;
	  abort();
  }

  //data_out.write_vtk (output);
  //data_out.outputFileType(output);

  //create pvtu record
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0){
    std::vector<std::string> filenames;
    for (unsigned int i=0;i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i) {
    	char vtuProcFileName[100];
    	sprintf(vtuProcFileName, "solution-%s.%u.%s", cycleAsString.str().c_str(),i,outputFileType);
    	filenames.push_back (vtuProcFileName);
    }
    std::ofstream master_output (pvtuFileName);

    data_out.write_pvtu_record (master_output, filenames);
  }
  pcout << "Output written to:" << pvtuFileName << "\n\n";

  //log time
  computing_timer.exit_section("matrixFreePDE: output"); 
}

#endif
