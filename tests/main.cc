#include <prismspf/utilities/utilities.h>

#include <catch2/catch_session.hpp>

int
main(int argc, char **argv)
{
  // NOTE: Always disallow threading
  prismspf::MPIInitFinalize mpi_init(argc, argv, 1);

  return Catch::Session().run(argc, argv);
}
