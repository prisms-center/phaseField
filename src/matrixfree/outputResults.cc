// outputResults() method for MatrixFreePDE class

#include <deal.II/numerics/data_out.h>

#include "../../include/matrixFreePDE.h"

// output results
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::outputResults()
{
  // log time
  computing_timer.enter_subsection("matrixFreePDE: output");

  // create DataOut object
  DataOut<dim> data_out;

  // loop over fields

  for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
    {
      // mark field as scalar/vector
      std::vector<DataComponentInterpretation::DataComponentInterpretation> dataType(
        fields[fieldIndex].numComponents,
        (fields[fieldIndex].type == SCALAR
           ? DataComponentInterpretation::component_is_scalar
           : DataComponentInterpretation::component_is_part_of_vector));
      // add field to data_out
      std::vector<std::string> solutionNames(fields[fieldIndex].numComponents,
                                             fields[fieldIndex].name.c_str());
      data_out.add_data_vector(*dofHandlersSet[fieldIndex],
                               *solutionSet[fieldIndex],
                               solutionNames,
                               dataType);
    }

  // Test section for outputting postprocessed fields
  // Currently there are hacks in place, using the matrixFreeObject, invM,
  // constraints, and DoFHandler as the primary variables Currently only works
  // for scalar fields. Require separate loop for vector-type post process
  // fields
  if (userInputs.postProcessingRequired)
    {
      std::vector<vectorType *> postProcessedSet;
      computePostProcessedFields(postProcessedSet);
#if (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR < 4)
      unsigned int invM_size = invMscalar.local_size();
      for (auto &field : postProcessedSet)
        {
          for (unsigned int dof = 0; dof < field->local_size(); ++dof)
            {
#else
      unsigned int invM_size = invMscalar.locally_owned_size();
      for (auto &field : postProcessedSet)
        {
          for (unsigned int dof = 0; dof < field->locally_owned_size(); ++dof)
            {
#endif
              field->local_element(dof) =
                invMscalar.local_element(dof % invM_size) * field->local_element(dof);
            }
          constraintsOtherSet[0]->distribute(*field);
          field->update_ghost_values();
        }

      // Integrate over selected post-processed fields and output them to the
      // screen and a text file
      std::ofstream output_file;

      if (userInputs.num_integrated_fields > 0)
        {
          if (first_integrated_var_output_complete)
            {
              output_file.open("integratedFields.txt", std::ios::app);
            }
          else
            {
              output_file.open("integratedFields.txt", std::ios::out);
            }
          output_file.precision(10);

          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            {
              output_file << currentTime;
            }

          for (const auto &[pp_index, pp_variable] : var_attributes.pp_attributes)
            {
              if (pp_variable.calc_integral)
                {
                  double integrated_field;
                  computeIntegral(integrated_field, pp_index, postProcessedSet);
                  pcout << "Integrated value of " << pp_variable.name << ": "
                        << integrated_field << std::endl;
                  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
                    {
                      output_file << "\t" << pp_variable.name << "\t" << integrated_field;
                    }
                  integrated_postprocessed_fields.at(
                    userInputs.integrated_field_indices[pp_index]) = integrated_field;
                }
            }
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            {
              output_file << std::endl;
            }
          output_file.close();
          first_integrated_var_output_complete = true;
        }

      // Add the postprocessed fields to data_out
      for (const auto &[fieldIndex, pp_variable] : var_attributes.pp_attributes)
        {
          // mark field as scalar/vector
          unsigned int components;
          if (userInputs.pp_varInfoList[fieldIndex].is_scalar)
            {
              components = 1;
              std::vector<DataComponentInterpretation::DataComponentInterpretation>
                dataType(components, DataComponentInterpretation::component_is_scalar);
              std::vector<std::string> solutionNames(components,
                                                     pp_variable.name.c_str());
              // add field to data_out
              data_out.add_data_vector(*dofHandlersSet[0],
                                       *postProcessedSet[fieldIndex],
                                       solutionNames,
                                       dataType);
            }
          else
            {
              components = dim;
              std::vector<DataComponentInterpretation::DataComponentInterpretation>
                                       dataType(components,
                         DataComponentInterpretation::component_is_part_of_vector);
              std::vector<std::string> solutionNames(components,
                                                     pp_variable.name.c_str());
              // add field to data_out
              data_out.add_data_vector(*dofHandlersSet[0],
                                       *postProcessedSet[fieldIndex],
                                       solutionNames,
                                       dataType);
            }
        }
    }

  data_out.build_patches(degree);

  // Defining snprintf with no warnings because we don't can about truncation
#define snprintf_nowarn(...) (snprintf(__VA_ARGS__) < 0 ? abort() : (void) 0)

  // write to results file
  // file name
  std::ostringstream cycleAsString;
  cycleAsString << std::setw(std::floor(std::log10(userInputs.totalIncrements)) + 1)
                << std::setfill('0') << currentIncrement;
  char baseFileName[100], vtuFileName[100];
  snprintf_nowarn(baseFileName,
                  sizeof(baseFileName),
                  "%s-%s",
                  userInputs.output_file_name.c_str(),
                  cycleAsString.str().c_str());
  snprintf_nowarn(vtuFileName,
                  sizeof(vtuFileName),
                  "%s.%u.%s",
                  baseFileName,
                  Utilities::MPI::this_mpi_process(MPI_COMM_WORLD),
                  userInputs.output_file_type.c_str());

  // Write to file in either vtu or vtk format
  if (userInputs.output_file_type == "vtu")
    {
      // Set flags to output the time and cycle number as part of the vtu file
      dealii::DataOutBase::VtkFlags flags;
      flags.time                = currentTime;
      flags.cycle               = currentIncrement;
      flags.print_date_and_time = true;
      data_out.set_flags(flags);

      if (userInputs.output_vtu_per_process)
        {
          // Write the results to separate files for each process
          std::ofstream output(vtuFileName);
          data_out.write_vtu(output);

          // Create pvtu record that can be used to stitch together the results
          // from all the processes
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            {
              std::vector<std::string> filenames;
              for (unsigned int i = 0;
                   i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
                   ++i)
                {
                  char vtuProcFileName[100];
                  snprintf_nowarn(vtuProcFileName,
                                  sizeof(vtuProcFileName),
                                  "%s-%s.%u.%s",
                                  userInputs.output_file_name.c_str(),
                                  cycleAsString.str().c_str(),
                                  i,
                                  userInputs.output_file_type.c_str());
                  filenames.emplace_back(vtuProcFileName);
                }
              char pvtuFileName[100];
              snprintf_nowarn(pvtuFileName,
                              sizeof(pvtuFileName),
                              "%s.p%s",
                              baseFileName,
                              userInputs.output_file_type.c_str());
              std::ofstream master_output(pvtuFileName);

              data_out.write_pvtu_record(master_output, filenames);
              pcout << "Output written to:" << pvtuFileName << "\n\n";
            }
        }
      else
        {
          // Write the results to a file shared between all processes
          char svtuFileName[100];
          snprintf_nowarn(svtuFileName,
                          sizeof(svtuFileName),
                          "%s.%s",
                          baseFileName,
                          userInputs.output_file_type.c_str());
          data_out.write_vtu_in_parallel(svtuFileName, MPI_COMM_WORLD);
          pcout << "Output written to:" << svtuFileName << "\n\n";
        }
    }
  else if (userInputs.output_file_type == "vtk")
    {
      // Write the results to separate files for each process
      std::ofstream output(vtuFileName);
      data_out.write_vtk(output);
      pcout << "Output written to:" << vtuFileName << "\n\n";
    }
  else
    {
      std::cerr << "PRISMS-PF Error: The parameter 'outputFileType' must be "
                   "either \"vtu\" or \"vtk\""
                << std::endl;
      abort();
    }

  // log time
  computing_timer.leave_subsection("matrixFreePDE: output");
}

#include "../../include/matrixFreePDE_template_instantiations.h"
