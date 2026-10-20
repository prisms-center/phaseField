// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>

#include <prismspf/utilities/assert.h>
#include <prismspf/utilities/utilities.h>

#include <prismspf/config.h>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali-manager.h>
#  include <caliper/cali.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

MPIInitFinalize::MPIInitFinalize(int                          &argc,
                                 char                       **&argv,
                                 [[maybe_unused]] unsigned int _max_n_threads)
  : dealii::Utilities::MPI::MPI_InitFinalize(argc,
                                             argv,
                                             _max_n_threads == 0 ? max_n_threads
                                                                 : _max_n_threads)

{
  // Search for the output and caliper arguments to use. Can't use getopt here since if
  // will effect subsequent calls.
  std::string output_dir     = "outputs/";
  std::string caliper_config = "runtime-report,mem.highwatermark";

  for (int i = 1; i < argc - 1; ++i)
    {
      if (std::string(argv[i]) == "-o")
        {
          output_dir = argv[i + 1];
        }
      else if (std::string(argv[i]) == "-c")
        {
          caliper_config = argv[i + 1];
        }
    }

  // Set log file
  Logger::set_file(output_dir + "simulation.log");

  // Restrict deal.II console printing
  dealii::deallog.depth_console(0);

  // Set the custom failure handler for libassert
  libassert::set_failure_handler(failure_handler);

#ifdef PRISMS_PF_WITH_CALIPER
  // Set the Caliper configuration
  mgr.add(caliper_config.c_str());

  // Check for configuration errors
  if (mgr.error())
    {
      std::cerr << "Caliper error: " << mgr.error_msg() << "\n" << std::flush;
    }

  // Start configured performance measurements, if any
  mgr.start();
#endif
}

MPIInitFinalize::~MPIInitFinalize()
{
#ifdef PRISMS_PF_WITH_CALIPER
  // Close and flush Caliper config manager. Must be done before finalizing MPI
  mgr.flush();
#endif
  // The dealii::Utilities::MPI::MPI_InitFinalize is implicitly destructed
}

PRISMS_PF_END_NAMESPACE
