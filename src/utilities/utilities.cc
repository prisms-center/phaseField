
#include <deal.II/base/logstream.h>
#include <deal.II/base/mpi.h>

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
  // Restrict deal.II console printing
  dealii::deallog.depth_console(0);

#ifdef PRISMS_PF_WITH_CALIPER
  // Add some useful defaults for Caliper
  // TODO: Make this configurable at the command line
  mgr.add("runtime-report,mem.highwatermark");

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
