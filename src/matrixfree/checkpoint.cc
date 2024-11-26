#include <deal.II/base/exceptions.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include "../../include/matrixFreePDE.h"
#include <filesystem>

#ifdef DEAL_II_WITH_ZLIB
#  include <zlib.h>
#endif

// Save a checkpoint
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::save_checkpoint()
{
  computing_timer.enter_subsection("matrixFreePDE: save_checkpoint");
  unsigned int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (my_id == 0)
    {
      // if we have previously written a snapshot, then keep the last
      // snapshot in case this one fails to save. Note: static variables
      // will only be initialized once per model run.
      static bool previous_snapshot_exists = (userInputs.resume_from_checkpoint == true);

      if (previous_snapshot_exists)
        {
          move_file("restart.mesh", "restart.mesh.old");
          move_file("restart.mesh.info", "restart.mesh.info.old");
          move_file("restart.time.info", "restart.time.info.old");
        }
      // from now on, we know that if we get into this
      // function again that a snapshot has previously
      // been written
      previous_snapshot_exists = true;
    }

  // save Triangulation and Solution vectors:
  {
    // Serializing all of the scalars together and all of the vectors together

    // First, get lists of scalar and vector fields
    std::vector<unsigned int> scalar_var_indices;
    std::vector<unsigned int> vector_var_indices;
    for (const auto &[index, variable] : var_attributes.attributes)
      {
        if (variable.var_type == SCALAR)
          {
            scalar_var_indices.push_back(index);
          }
        else
          {
            vector_var_indices.push_back(index);
          }
      }

    // Second, build one solution set list for scalars and one for vectors
    std::vector<const vectorType *> solSet_transfer_scalars;
    std::vector<const vectorType *> solSet_transfer_vectors;
    for (const auto &[index, variable] : var_attributes.attributes)
      {
        if (variable.var_type == SCALAR)
          {
            solSet_transfer_scalars.push_back(solutionSet[index]);
          }
        else
          {
            solSet_transfer_vectors.push_back(solutionSet[index]);
          }
      }

    // Finally, save the triangulation and the solutionTransfer objects
    if (scalar_var_indices.size() > 0 && vector_var_indices.size() == 0)
      {
        parallel::distributed::SolutionTransfer<dim, vectorType> system_trans_scalars(
          *dofHandlersSet[scalar_var_indices[0]]);
        system_trans_scalars.prepare_for_serialization(solSet_transfer_scalars);

        triangulation.save("restart.mesh");
      }
    else if (scalar_var_indices.size() == 0 && vector_var_indices.size() > 0)
      {
        parallel::distributed::SolutionTransfer<dim, vectorType> system_trans_vectors(
          *dofHandlersSet[vector_var_indices[0]]);
        system_trans_vectors.prepare_for_serialization(solSet_transfer_vectors);

        triangulation.save("restart.mesh");
      }
    else
      {
        parallel::distributed::SolutionTransfer<dim, vectorType> system_trans_scalars(
          *dofHandlersSet[scalar_var_indices[0]]);
        system_trans_scalars.prepare_for_serialization(solSet_transfer_scalars);

        parallel::distributed::SolutionTransfer<dim, vectorType> system_trans_vectors(
          *dofHandlersSet[vector_var_indices[0]]);
        system_trans_vectors.prepare_for_serialization(solSet_transfer_vectors);

        triangulation.save("restart.mesh");
      }
  }

  // Save information about the current increment and current time
  if (my_id == 0)
    {
      std::ofstream time_info_file;
      time_info_file.open("restart.time.info");
      time_info_file << currentIncrement << " (currentIncrement)\n";
      time_info_file << currentTime << " (currentTime)\n";
      time_info_file.close();
    }

  pcout << "*** Checkpoint created! ***\n\n";
  computing_timer.leave_subsection("matrixFreePDE: save_checkpoint");
}

// Load from a previously created checkpoint
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::load_checkpoint_triangulation()
{
  // First check existence of the two restart files for the mesh and field
  // variables
  verify_checkpoint_file_exists("restart.mesh");
  verify_checkpoint_file_exists("restart.mesh.info");

  pcout << "\n*** Resuming from a checkpoint! ***\n\n";

  try
    {
      triangulation.load("restart.mesh");
    }
  catch (...)
    {
      AssertThrow(false,
                  ExcMessage("PRISMS-PF Error: Cannot open snapshot mesh file "
                             "or read the triangulation stored there."));
    }
}

// Load from a previously saved checkpoint
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::load_checkpoint_fields()
{
  // Serializing all of the scalars together and all of the vectors together

  // First, get lists of scalar and vector fields
  std::vector<unsigned int> scalar_var_indices;
  std::vector<unsigned int> vector_var_indices;
  for (const auto &[index, variable] : var_attributes.attributes)
    {
      if (variable.var_type == SCALAR)
        {
          scalar_var_indices.push_back(index);
        }
      else
        {
          vector_var_indices.push_back(index);
        }
    }

  // Second, build one solution set list for scalars and one for vectors
  std::vector<vectorType *> solSet_transfer_scalars;
  std::vector<vectorType *> solSet_transfer_vectors;
  for (const auto &[index, variable] : var_attributes.attributes)
    {
      if (variable.var_type == SCALAR)
        {
          solSet_transfer_scalars.push_back(solutionSet[index]);
        }
      else
        {
          solSet_transfer_vectors.push_back(solutionSet[index]);
        }
    }

  // Finally, deserialize the fields to the solSet_transfer objects, which
  // contain pointers to solutionSet
  if (scalar_var_indices.size() > 0)
    {
      parallel::distributed::SolutionTransfer<dim, vectorType> system_trans_scalars(
        *dofHandlersSet[scalar_var_indices[0]]);
      system_trans_scalars.deserialize(solSet_transfer_scalars);
    }
  if (vector_var_indices.size() > 0)
    {
      parallel::distributed::SolutionTransfer<dim, vectorType> system_trans_vectors(
        *dofHandlersSet[vector_var_indices[0]]);
      system_trans_vectors.deserialize(solSet_transfer_vectors);
    }
}

// Load from a previously saved checkpoint
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::load_checkpoint_time_info()
{
  // Make sure that restart.time.info exists
  verify_checkpoint_file_exists("restart.time.info");

  std::ifstream time_info_file;
  time_info_file.open("restart.time.info");
  std::string line;
  std::getline(time_info_file, line);
  line.erase(line.end() - 19, line.end());
  currentIncrement = dealii::Utilities::string_to_int(line);

  std::getline(time_info_file, line);
  line.erase(line.end() - 14, line.end());
  currentTime = dealii::Utilities::string_to_double(line);
  time_info_file.close();
}

// Move/rename a checkpoint file
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::move_file(const std::string &old_name,
                                      const std::string &new_name)
{
  try
    {
      std::filesystem::rename(old_name, new_name);
    }
  catch (const std::filesystem::filesystem_error &error)
    {
      if (std::filesystem::exists(new_name))
        {
          bool is_removed = std::filesystem::remove(new_name);
          AssertThrow(is_removed,
                      dealii::ExcMessage(
                        "Unable to remove file: " + new_name +
                        ", although it seems to exist. The error code is " +
                        error.what()));
        }

      try
        {
          std::filesystem::rename(old_name, new_name);
        }
      catch (const std::filesystem::filesystem_error &rename_error)
        {
          AssertThrow(false,
                      ExcMessage("Unable to rename file: " + old_name + " -> " +
                                 new_name + ". Error: " + rename_error.what()));
        }
    }
}

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::verify_checkpoint_file_exists(const std::string &filename)
{
  std::ifstream in(filename);
  if (!in)
    {
      AssertThrow(false,
                  ExcMessage(std::string("PRISMS-PF Error: You are trying to "
                                         "restart a previous computation, "
                                         "but the restart file <") +
                             filename + "> does not appear to exist!"));
    }
}
