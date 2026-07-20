#
# Find MPI and run some checks
#

find_package(MPI REQUIRED COMPONENTS CXX)

if(NOT MPI_CXX_FOUND)
  message(
    FATAL_ERROR
    "MPI CXX component not found. Please make sure your MPI installation includes C++ support."
  )
endif()
if(NOT MPI_CXX_COMPILER)
  message(
    FATAL_ERROR
    "MPI CXX compiler wrapper not found. Please make sure you have `mpicxx`"
  )
endif()

# Find the MPI CXX version if not already populated
if(NOT MPI_CXX_VERSION)
  file(READ "${MPI_CXX_INCLUDE_DIRS}/mpi.h" _mpi_header_content)

  string(REGEX MATCH "#define MPI_MAJOR_VERSION ([0-9]+)" _ "${_mpi_header_content}")
  set(_mpi_major ${CMAKE_MATCH_1})

  string(REGEX MATCH "#define MPI_MINOR_VERSION ([0-9]+)" _ "${_mpi_header_content}")
  set(_mpi_major ${CMAKE_MATCH_1})

  set(MPI_CXX_VERSION "${_mpi_major}.${_mpi_minor}")
endif()

# Add MPI to the Release and Debug lists
prisms_pf_add_dependency_target(MPI::MPI_CXX DEBUG PUBLIC)
prisms_pf_add_dependency_target(MPI::MPI_CXX RELEASE PUBLIC)
