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
      for (unsigned int fieldIndex = 0; fieldIndex < postProcessedSet.size();
           fieldIndex++)
        {
          for (unsigned int dof = 0; dof < postProcessedSet[fieldIndex]->local_size();
               ++dof)
            {
#else
      unsigned int invM_size = invMscalar.locally_owned_size();
      for (unsigned int fieldIndex = 0; fieldIndex < postProcessedSet.size();
           fieldIndex++)
        {
          for (unsigned int dof = 0;
               dof < postProcessedSet[fieldIndex]->locally_owned_size();
               ++dof)
            {
#endif
              postProcessedSet[fieldIndex]->local_element(dof) =
                invMscalar.local_element(dof % invM_size) *
                postProcessedSet[fieldIndex]->local_element(dof);
            }
          constraintsOtherSet[0]->distribute(*postProcessedSet[fieldIndex]);
          postProcessedSet[fieldIndex]->update_ghost_values();
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

          for (unsigned int i = 0; i < userInputs.pp_number_of_variables; i++)
            {
              if (userInputs.pp_calc_integral[i])
                {
                  double integrated_field;
                  computeIntegral(integrated_field, i, postProcessedSet);
                  pcout << "Integrated value of "
                        << userInputs.pp_var_name[userInputs.integrated_field_indices[i]]
                        << ": " << integrated_field << std::endl;
                  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
                    {
                      output_file
                        << "\t"
                        << userInputs.pp_var_name[userInputs.integrated_field_indices[i]]
                        << "\t" << integrated_field;
                    }
                  integrated_postprocessed_fields.at(
                    userInputs.integrated_field_indices[i]) = integrated_field;
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
      for (unsigned int fieldIndex = 0; fieldIndex < userInputs.pp_number_of_variables;
           fieldIndex++)
        {
          // mark field as scalar/vector
          unsigned int components;
          if (userInputs.pp_varInfoList[fieldIndex].is_scalar)
            {
              components = 1;
              std::vector<DataComponentInterpretation::DataComponentInterpretation>
                dataType(components, DataComponentInterpretation::component_is_scalar);
              std::vector<std::string> solutionNames(
                components,
                userInputs.pp_var_name[fieldIndex].c_str());
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
              std::vector<std::string> solutionNames(
                components,
                userInputs.pp_var_name[fieldIndex].c_str());
              // add field to data_out
              data_out.add_data_vector(*dofHandlersSet[0],
                                       *postProcessedSet[fieldIndex],
                                       solutionNames,
                                       dataType);
            }
        }
    }

  data_out.build_patches(degree);

  // write to results file
  // file name
  std::ostringstream cycleAsString;
  cycleAsString << std::setw(std::floor(std::log10(userInputs.totalIncrements)) + 1)
                << std::setfill('0') << currentIncrement;
  char baseFileName[100], vtuFileName[100];
  snprintf(baseFileName,
           sizeof(baseFileName),
           "%s-%s",
           userInputs.output_file_name.c_str(),
           cycleAsString.str().c_str());
  snprintf(vtuFileName,
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
                  snprintf(vtuProcFileName,
                           sizeof(vtuProcFileName),
                           "%s-%s.%u.%s",
                           userInputs.output_file_name.c_str(),
                           cycleAsString.str().c_str(),
                           i,
                           userInputs.output_file_type.c_str());
                  filenames.push_back(vtuProcFileName);
                }
              char pvtuFileName[100];
              snprintf(pvtuFileName,
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
          snprintf(svtuFileName,
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
