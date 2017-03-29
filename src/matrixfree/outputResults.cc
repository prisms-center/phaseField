//outputResults() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

//output results
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::outputResults() const {
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
  
  // Test section for outputting postprocessed fields
  PostProcessor<dim,degree> post_processor(userInputs);
  std::vector<vectorType*> postProcessedSet;
  post_processor.computePostProcessedFields(matrixFreeObject,solutionSet,postProcessedSet);
  unsigned int invM_size = invM.local_size();
  for(unsigned int fieldIndex=0; fieldIndex<solutionSet.size(); fieldIndex++){
	  for (unsigned int dof=0; dof<solutionSet[fieldIndex]->local_size(); ++dof){
		  postProcessedSet[fieldIndex]->local_element(dof)=			\
				  invM.local_element(dof%invM_size)*postProcessedSet[fieldIndex]->local_element(dof);
	  }
  }

  std::vector<std::string> pp_fields;
  pp_fields.push_back("n_pp");
  for(unsigned int fieldIndex=0; fieldIndex<1; fieldIndex++){
	  //mark field as scalar/vector
	  std::vector<DataComponentInterpretation::DataComponentInterpretation> dataType \
	  (fields[fieldIndex].numComponents,				\
			  (fields[fieldIndex].type==SCALAR ?				\
					  DataComponentInterpretation::component_is_scalar:		\
					  DataComponentInterpretation::component_is_part_of_vector));
	  //add field to data_out
	  std::vector<std::string> solutionNames (fields[fieldIndex].numComponents, pp_fields[fieldIndex].c_str());
	  data_out.add_data_vector(*dofHandlersSet[fieldIndex], *postProcessedSet[fieldIndex], solutionNames, dataType);
  }

  data_out.build_patches (degree);
  
  //write to results file
  //file name
  std::ostringstream cycleAsString;
  cycleAsString << std::setw(std::ceil(std::log10(userInputs.totalIncrements))+1) << std::setfill('0') << currentIncrement;
  char vtuFileName[100], pvtuFileName[100];
  sprintf(vtuFileName, "solution-%s.%u.%s", cycleAsString.str().c_str(),Utilities::MPI::this_mpi_process(MPI_COMM_WORLD),userInputs.output_file_type.c_str());
  sprintf(pvtuFileName, "solution-%s.p%s", cycleAsString.str().c_str(),userInputs.output_file_type.c_str());
  std::ofstream output (vtuFileName);

  //write to file
  if (userInputs.output_file_type == "vtu"){
	  data_out.write_vtu (output);
  }
  else if (userInputs.output_file_type == "vtk"){
	  data_out.write_vtk (output);
  }
  else {
	  std::cout << "PRISMS-PF Error: The parameter 'outputFileType' must be either \"vtu\" or \"vtk\"" << std::endl;
	  abort();
  }

  //create pvtu record
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0){
    std::vector<std::string> filenames;
    for (unsigned int i=0;i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i) {
    	char vtuProcFileName[100];
    	sprintf(vtuProcFileName, "solution-%s.%u.%s", cycleAsString.str().c_str(),i,userInputs.output_file_type.c_str());
    	filenames.push_back (vtuProcFileName);
    }
    std::ofstream master_output (pvtuFileName);

    data_out.write_pvtu_record (master_output, filenames);
  }
  pcout << "Output written to:" << pvtuFileName << "\n\n";

  //log time
  computing_timer.exit_section("matrixFreePDE: output"); 
}

#include "../../include/matrixFreePDE_template_instantiations.h"
