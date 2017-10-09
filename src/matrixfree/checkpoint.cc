#include "../../include/matrixFreePDE.h"

#ifdef DEAL_II_WITH_ZLIB
#  include <zlib.h>
#endif

// Save a checkpoint
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::save_checkpoint(){

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
            //move_file ("restart.resume.z","restart.resume.z.old");
        }
        // from now on, we know that if we get into this
        // function again that a snapshot has previously
        // been written
        previous_snapshot_exists = true;
    }

    // save Triangulation and Solution vectors:
    {
    //   std::vector<const LinearAlgebra::BlockVector *> x_system (3);
    //   x_system[0] = &solution;
    //   x_system[1] = &old_solution;
    //   x_system[2] = &old_old_solution;
      //

    std::vector<vectorType*> solutionSet_for_transfer = solutionSet;

    parallel::distributed::SolutionTransfer<dim, vectorType> system_trans (*dofHandlersSet[0]);
    system_trans.prepare_serialization (*solutionSet_for_transfer[0]);

      //
    //   system_trans.prepare_serialization (x_system);

      //signals.pre_checkpoint_store_user_data(triangulation);

      triangulation.save ("restart.mesh");
    }

    // save general information This calls the serialization functions on all
    // processes (so that they can take additional action, if necessary, see
    // the manual) but only writes to the restart file on process 0
    {
      std::ostringstream oss;

      // serialize into a stringstream
      boost::archive::binary_oarchive oa (oss);
      oa << (*this);

      // compress with zlib and write to file on the root processor
#ifdef DEAL_II_WITH_ZLIB
      if (my_id == 0)
        {
          uLongf compressed_data_length = compressBound (oss.str().length());
          std::vector<char *> compressed_data (compressed_data_length);
          int err = compress2 ((Bytef *) &compressed_data[0],
                               &compressed_data_length,
                               (const Bytef *) oss.str().data(),
                               oss.str().length(),
                               Z_BEST_COMPRESSION);
          (void)err;
          Assert (err == Z_OK, ExcInternalError());

          // build compression header
          const uint32_t compression_header[4]
            = { 1,                                   /* number of blocks */
                (uint32_t)oss.str().length(), /* size of block */
                (uint32_t)oss.str().length(), /* size of last block */
                (uint32_t)compressed_data_length
              }; /* list of compressed sizes of blocks */

          std::ofstream f ("restart.resume.z");
          f.write((const char *)compression_header, 4 * sizeof(compression_header[0]));
          f.write((char *)&compressed_data[0], compressed_data_length);
        }
#else
      AssertThrow (false,
                   ExcMessage ("You need to have deal.II configured with the `libz' "
                               "option to support checkpoint/restart, but deal.II "
                               "did not detect its presence when you called `cmake'."));
#endif

     }

    pcout << "*** Snapshot created!" << std::endl << std::endl;

}


// Load from a previously saved checkpoint
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::load_checkpoint_triangulation(){

    // first check existence of the two restart files
    {
      const std::string filename = "restart.mesh";
      std::ifstream in (filename.c_str());
      if (!in)
        AssertThrow (false,
                     ExcMessage (std::string("You are trying to restart a previous computation, "
                                             "but the restart file <")
                                 +
                                 filename
                                 +
                                 "> does not appear to exist!"));
    }
    /*{
      const std::string filename = parameters.output_directory + "restart.resume.z";
      std::ifstream in (filename.c_str());
      if (!in)
        AssertThrow (false,
                     ExcMessage (std::string("You are trying to restart a previous computation, "
                                             "but the restart file <")
                                 +
                                 filename
                                 +
                                 "> does not appear to exist!"));
    }*/

    pcout << "*** Resuming from snapshot!" << std::endl << std::endl;

    try
    {
        triangulation.load ("restart.mesh");
    }
    catch (...)
    {
        AssertThrow(false, ExcMessage("Cannot open snapshot mesh file or read the triangulation stored there."));
    }

}

// Load from a previously saved checkpoint
template <int dim, int degree>
void MatrixFreePDE<dim,degree>::load_checkpoint_fields(){

    parallel::distributed::SolutionTransfer<dim, vectorType> system_trans (*dofHandlersSet[0]);
    system_trans.deserialize (*solutionSet[0]);

    // // first check existence of the two restart files
    // {
    //   const std::string filename = "restart.mesh";
    //   std::ifstream in (filename.c_str());
    //   if (!in)
    //     AssertThrow (false,
    //                  ExcMessage (std::string("You are trying to restart a previous computation, "
    //                                          "but the restart file <")
    //                              +
    //                              filename
    //                              +
    //                              "> does not appear to exist!"));
    // }
    // /*{
    //   const std::string filename = parameters.output_directory + "restart.resume.z";
    //   std::ifstream in (filename.c_str());
    //   if (!in)
    //     AssertThrow (false,
    //                  ExcMessage (std::string("You are trying to restart a previous computation, "
    //                                          "but the restart file <")
    //                              +
    //                              filename
    //                              +
    //                              "> does not appear to exist!"));
    // }*/
    //
    // pcout << "*** Resuming from snapshot!" << std::endl << std::endl;
    //
    // try
    // {
    //     triangulation.load ("restart.mesh");
    // }
    // catch (...)
    // {
    //     AssertThrow(false, ExcMessage("Cannot open snapshot mesh file or read the triangulation stored there."));
    // }


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
