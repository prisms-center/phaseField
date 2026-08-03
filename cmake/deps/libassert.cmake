#
# Find libassert and run some checks
#

# With libassert we always want to vendor it

# Add the external projects
prisms_pf_add_external_project(
  libassert
  https://github.com/jeremy-rifkin/libassert.git
  v2.2.1
  libassert
  libcpptrace
  libzstd
  libdwarf
)

# Create the libraries
prisms_pf_add_external_library(libassert libassert)

# Add transitive dependencies and link them to libassert
# NOTE: We have to do this to make use of libassert when linking
# directly to the prisms_pf target (e.g., unit tests).
find_package(Threads)
# dwarf
prisms_pf_import_archive(imported_dwarf "${CMAKE_BINARY_DIR}/_deps/libassert" libdwarf)
# zstd
prisms_pf_import_archive(imported_zstd "${CMAKE_BINARY_DIR}/_deps/libassert" libzstd)
# cpptrace
prisms_pf_import_archive(
  imported_cpptrace
  "${CMAKE_BINARY_DIR}/_deps/libassert"
  libcpptrace
)
set_target_properties(
  imported_cpptrace
  PROPERTIES
    INTERFACE_LINK_LIBRARIES
      "imported_dwarf;imported_zstd;Threads::Threads;${CMAKE_DL_LIBS}"
)
if(CMAKE_BUILD_TYPE STREQUAL "DebugRelease")
  set_target_properties(
    imported_cpptrace_debug
    PROPERTIES
      INTERFACE_LINK_LIBRARIES
        "imported_dwarf_debug;imported_zstd_debug;Threads::Threads;${CMAKE_DL_LIBS}"
  )
endif()
set_target_properties(
  imported_libassert
  PROPERTIES
    INTERFACE_LINK_LIBRARIES
      imported_cpptrace
)
if(CMAKE_BUILD_TYPE STREQUAL "DebugRelease")
  set_target_properties(
    imported_libassert_debug
    PROPERTIES
      INTERFACE_LINK_LIBRARIES
        imported_cpptrace_debug
  )
endif()

# Create install rules
prisms_pf_install_external_library(libassert)

# Add libassert to the Release and Debug lists
if(CMAKE_BUILD_TYPE STREQUAL "DebugRelease")
  prisms_pf_add_dependency_target(imported_libassert_debug DEBUG PUBLIC)
  prisms_pf_add_dependency_target(imported_libassert RELEASE PUBLIC)
else()
  prisms_pf_add_dependency_target(imported_libassert "${CMAKE_BUILD_TYPE}" PUBLIC)
endif()
