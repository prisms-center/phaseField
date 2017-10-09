#include "../../include/matrixFreePDE.h"

#ifdef DEAL_II_WITH_ZLIB
#  include <zlib.h>
#endif

// Save a checkpoint
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::save_checkpoint(){
    computing_timer.enter_section("matrixFreePDE: save_checkpoint");
    unsigned int my_id = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

    if (my_id == 0)
    {
        // if we have previously written a snapshot, then keep the last
        // snapshot in case this one fails to save. Note: static variables
        // will only be initialized once per model run.
        static bool previous_snapshot_exists = (userInputs.resume_from_checkpoint == true);

        if (previous_snapshot_exists == true)
        {
            move_file ("restart.mesh","restart.mesh.old");
            move_file ("restart.mesh.info","restart.mesh.info.old");
            move_file ("restart.time.info","restart.time.info.old");
        }
        // from now on, we know that if we get into this
        // function again that a snapshot has previously
        // been written
        previous_snapshot_exists = true;
    }

    // save Triangulation and Solution vectors:
    {

    std::vector<vectorType*> solutionSet_for_transfer = solutionSet;

    parallel::distributed::SolutionTransfer<dim, vectorType> system_trans (*dofHandlersSet[0]);
    system_trans.prepare_serialization (*solutionSet_for_transfer[0]);

      triangulation.save ("restart.mesh");
    }

    // Save information about the current increment and current time
    if (my_id == 0){
        std::ofstream time_info_file;
        time_info_file.open("restart.time.info");
        time_info_file << currentIncrement << " (currentIncrement)\n";
        time_info_file << currentTime << " (currentTime)\n";
        time_info_file.close();
    }

    pcout << "*** Checkpoint created!" << std::endl << std::endl;
    computing_timer.exit_section("matrixFreePDE: save_checkpoint");

}


// Load from a previously created checkpoint
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::load_checkpoint_triangulation(){

    // first check existence of the two restart files
    {
      const std::string filename = "restart.mesh";
      std::ifstream in (filename.c_str());
      if (!in)
        AssertThrow (false,
                     ExcMessage (std::string("PRISMS-PF Error: You are trying to restart a previous computation, "
                                             "but the restart file <")
                                 +
                                 filename
                                 +
                                 "> does not appear to exist!"));
    }

    pcout << "*** Resuming from a checkpoint!" << std::endl << std::endl;

    try
    {
        triangulation.load ("restart.mesh");
    }
    catch (...)
    {
        AssertThrow(false, ExcMessage("PRISMS-PF Error: Cannot open snapshot mesh file or read the triangulation stored there."));
    }

}

// Load from a previously saved checkpoint
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::load_checkpoint_fields(){

    parallel::distributed::SolutionTransfer<dim, vectorType> system_trans (*dofHandlersSet[0]);
    system_trans.deserialize (*solutionSet[0]);

}

// Load from a previously saved checkpoint
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::load_checkpoint_time_info(){

    // Make sure that restart.time.info exists
    {
        const std::string filename = "restart.time.info";

        std::ifstream in (filename);
        if (!in){
            AssertThrow (false,
                ExcMessage (std::string("PRISMS-PF Error: You are trying to restart a previous computation, "
                "but the restart file <")
                +
                filename
                +
                "> does not appear to exist!"));
            }
    }

    std::ifstream time_info_file;
    time_info_file.open("restart.time.info");
    std::string line;
    std::getline(time_info_file, line);
    line.erase(line.end()-19,line.end());
    currentIncrement = dealii::Utilities::string_to_int(line);

    std::getline(time_info_file, line);
    line.erase(line.end()-14,line.end());
    currentTime = dealii::Utilities::string_to_double(line);
    time_info_file.close();

}


// Move/rename a checkpoint file
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::move_file (const std::string &old_name, const std::string &new_name){

    int error = system (("mv " + old_name + " " + new_name).c_str());

    // If the above call failed, e.g. because there is no command-line
    // available, try with internal functions.
    if (error != 0)
    {
        std::ifstream ifile(new_name);
        if (static_cast<bool>(ifile))
        {
            error = remove(new_name.c_str());
            AssertThrow (error == 0, ExcMessage(std::string ("Unable to remove file: "
            + new_name
            + ", although it seems to exist. "
            + "The error code is "
            + dealii::Utilities::to_string(error) + ".")));
        }

        error = rename(old_name.c_str(),new_name.c_str());
        AssertThrow (error == 0, ExcMessage(std::string ("Unable to rename files: ")
        +
        old_name + " -> " + new_name
        + ". The error code is "
        + dealii::Utilities::to_string(error) + "."));
    }
}

#include "../../include/matrixFreePDE_template_instantiations.h"
