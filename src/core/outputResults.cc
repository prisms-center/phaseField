#include <deal.II/numerics/data_out.h>

#include <cmath>
#include <core/matrixFreePDE.h>
#include <fstream>
#include <sstream>

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
      std::vector<dealii::LinearAlgebra::distributed::Vector<double> *> postProcessedSet;
      computePostProcessedFields(postProcessedSet);
      unsigned int invM_size = invMscalar.locally_owned_size();
      for (auto &field : postProcessedSet)
        {
          for (unsigned int dof = 0; dof < field->locally_owned_size(); ++dof)
            {
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

          for (const auto &[pp_index, pp_variable] : pp_attributes)
            {
              if (pp_variable.calc_integral)
                {
                  double integrated_field = NAN;
                  computeIntegral(integrated_field, pp_index, postProcessedSet);
                  pcout << "Integrated value of " << pp_variable.name << ": "
                        << integrated_field << std::endl;
                  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
                    {
                      output_file << "\t" << pp_variable.name << "\t" << integrated_field;
                    }
                }
            }
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            {
              output_file << "\n";
            }
          output_file.close();
          first_integrated_var_output_complete = true;
        }

      // Add the postprocessed fields to data_out
      for (const auto &[fieldIndex, pp_variable] : pp_attributes)
        {
          // mark field as scalar/vector
          unsigned int components = 0;
          if (userInputs.pp_varInfoList[fieldIndex].is_scalar)
            {
              components = 1;
              std::vector<DataComponentInterpretation::DataComponentInterpretation>
                dataType(components, DataComponentInterpretation::component_is_scalar);
              std::vector<std::string> solutionNames(components, pp_variable.name);
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
              std::vector<std::string> solutionNames(components, pp_variable.name);
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

  std::ostringstream baseFileNameStream;
  baseFileNameStream << userInputs.output_file_name << "-" << cycleAsString.str();
  std::string baseFileName = baseFileNameStream.str();

  std::ostringstream vtuFileNameStream;
  vtuFileNameStream << baseFileName << "."
                    << Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) << "."
                    << userInputs.output_file_type;
  std::string vtuFileName = vtuFileNameStream.str();

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
                  std::ostringstream vtuProcFileNameStream;
                  vtuProcFileNameStream << userInputs.output_file_name << "-"
                                        << cycleAsString.str() << "." << i << "."
                                        << userInputs.output_file_type;
                  std::string vtuProcFileName = vtuProcFileNameStream.str();

                  filenames.emplace_back(vtuProcFileName);
                }

              std::ostringstream pvtuFileNameStream;
              pvtuFileNameStream << baseFileName << ".p" << userInputs.output_file_type;
              std::string pvtuFileName = pvtuFileNameStream.str();

              std::ofstream master_output(pvtuFileName);

              data_out.write_pvtu_record(master_output, filenames);
              pcout << "Output written to:" << pvtuFileName << "\n\n";
            }
        }
      else
        {
          // Write the results to a file shared between all processes
          std::ostringstream svtuFileNameStream;
          svtuFileNameStream << baseFileName << "." << userInputs.output_file_type;
          std::string svtuFileName = svtuFileNameStream.str();

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
                   "either \"vtu\" or \"vtk\"\n";
      abort();
    }

  // log time
  computing_timer.leave_subsection("matrixFreePDE: output");
}

template class MatrixFreePDE<2, 1>;
template class MatrixFreePDE<3, 1>;

template class MatrixFreePDE<2, 2>;
template class MatrixFreePDE<3, 2>;

template class MatrixFreePDE<3, 3>;
template class MatrixFreePDE<2, 3>;

template class MatrixFreePDE<3, 4>;
template class MatrixFreePDE<2, 4>;

template class MatrixFreePDE<3, 5>;
template class MatrixFreePDE<2, 5>;

template class MatrixFreePDE<3, 6>;
template class MatrixFreePDE<2, 6>;