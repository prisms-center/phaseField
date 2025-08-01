#
# Try and find the deal.II library
#

set(DEAL_II_DIR "" CACHE PATH "An optional hint to a deal.II directory")
set_if_empty(DEAL_II_DIR "$ENV{DEAL_II_DIR}")

find_package(deal.II 9.6.0 QUIET HINTS ${DEAL_II_DIR})
if(NOT ${deal.II_FOUND})
    message(
        FATAL_ERROR
        "\n*** Could not find a recent version of deal.II. ***\n"
        "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake "
        "or set an environment variable \"DEAL_II_DIR\" that contains a path to a "
        "recent version of deal.II."
    )
endif()
message(
    STATUS
    "Found deal.II version ${DEAL_II_PACKAGE_VERSION} at '${deal.II_DIR}'"
)

# Check the deal.II build types match the ones we are requesting
string(FIND "${DEAL_II_BUILD_TYPE}" "${CMAKE_BUILD_TYPE}" _pos)
if(_pos EQUAL -1)
    message(
        FATAL_ERROR
        "Mismatch between the build types that deal.II was built with and "
        "the ones being requested. PRISMS-PF is tryng to build with ${CMAKE_BUILD_TYPE}, but "
        "deal.II was built with ${DEAL_II_BUILD_TYPE}."
    )
endif()

# Required deal.II features
set(DEAL_II_INSTALL_VALID ON)

if(NOT DEAL_II_WITH_LAPACK)
    message(SEND_ERROR "\n**deal.II was built without support for lapack!\n")
    set(DEAL_II_INSTALL_VALID OFF)
endif()

if(NOT DEAL_II_WITH_P4EST)
    message(SEND_ERROR "\n**deal.II was built without support for p4est!\n")
    set(DEAL_II_INSTALL_VALID OFF)
endif()

if(NOT DEAL_II_WITH_MPI)
    message(SEND_ERROR "\n**deal.II was built without support for MPI!\n")
    set(DEAL_II_INSTALL_VALID OFF)
endif()

if(NOT DEAL_II_INSTALL_VALID)
    message(
        FATAL_ERROR
        "\nPRISMS-PF requires a deal.II installation with certain features enabled!\n"
    )
endif()

# Optional deal.II features

# zlib should be required by p4est in deal.II, so its on by default
if(PRISMS_PF_WITH_ZLIB)
    if(DEAL_II_WITH_ZLIB)
        message(STATUS "  Found deal.II installation with zlib")
    else()
        message(
            FATAL_ERROR
            "deal.II installation with zlib not found. Disable PRISMS_PF_WITH_ZLIB or recompile deal.II with zlib."
        )
    endif()
endif()

if(PRISMS_PF_AUTODETECTION AND DEAL_II_WITH_HDF5)
    set(PRISMS_PF_WITH_HDF5 ON)
endif()
if(PRISMS_PF_WITH_HDF5)
    if(DEAL_II_WITH_HDF5)
        message(STATUS "  Found deal.II installation with HDF5")
    else()
        message(
            FATAL_ERROR
            "deal.II installation with HDF5 not found. Disable PRISMS_PF_WITH_HDF5 or recompile deal.II with HDF5."
        )
    endif()
endif()

if(PRISMS_PF_AUTODETECTION AND DEAL_II_WITH_SUNDIALS)
    set(PRISMS_PF_WITH_SUNDIALS ON)
endif()
if(PRISMS_PF_WITH_SUNDIALS)
    if(DEAL_II_WITH_SUNDIALS)
        message(STATUS "  Found deal.II installation with SUNDIALS")
    else()
        message(
            FATAL_ERROR
            "deal.II installation with SUNDIALS not found. Disable PRISMS_PF_WITH_SUNDIALS or recompile deal.II with SUNDIALS."
        )
    endif()
endif()

# To enable CUDA in PRISMS-PF, deal.II must be built with Kokkos
set(KOKKOS_CUDA_BACKEND
    OFF
    CACHE BOOL
    "Whether the installed Kokkos version has a CUDA backend."
)
mark_as_advanced(KOKKOS_CUDA_BACKEND)
if(DEAL_II_WITH_KOKKOS)
    find_package(Kokkos QUIET HINTS ${KOKKOS_DIR} $ENV{KOKKOS_DIR})
    if(Kokkos_FOUND AND Kokkos_DEVICES MATCHES "CUDA")
        set(KOKKOS_CUDA_BACKEND ON)
    endif()
endif()
if(PRISMS_PF_AUTODETECTION AND KOKKOS_CUDA_BACKEND)
    set(PRISMS_PF_WITH_CUDA ON)
endif()
if(PRISMS_PF_WITH_CUDA)
    if(KOKKOS_CUDA_BACKEND)
        message(STATUS "  Found Kokkos installation with CUDA")
    else()
        message(
            FATAL_ERROR
            "Kokkos installation with CUDA not found. Disable PRISMS_PF_WITH_CUDA or recompile deal.II with Kokkos that has a CUDA backend."
        )
    endif()
endif()

# Grab some compiler flags that would otherwise be provided by deal_ii_initialize_cached_variables()
message(STATUS ${DEAL_II_CXX_COMPILER})
message(STATUS ${DEAL_II_CXX_FLAGS})
message(STATUS ${DEAL_II_CXX_FLAGS_DEBUG})
message(STATUS ${DEAL_II_CXX_FLAGS_RELEASE})

message(STATUS ${DEAL_II_LINKER_FLAGS})
message(STATUS ${DEAL_II_LINKER_FLAGS_DEBUG})
message(STATUS ${DEAL_II_LINKER_FLAGS_RELEASE})

message(STATUS ${DEAL_II_WARNING_FLAGS})
