#
# Find libassert and run some checks
#

# With libassert we always want to vendor it

function(add_external_project NAME GIT_REPO GIT_TAG)
  if(NOT NAME IN_LIST PRISMS_PF_VENDORED_PACKAGES)
    message(FATAL_ERROR "Invalid vendored package name")
    return()
  endif()

  ExternalProject_Add(
    ${NAME}
    GIT_REPOSITORY ${GIT_REPO}
    GIT_TAG ${GIT_TAG}
    PREFIX "${CMAKE_BINARY_DIR}/_deps/${NAME}"
    INSTALL_DIR "${CMAKE_BINARY_DIR}/_deps/${NAME}/install"
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ${ARGN}
  )
  ExternalProject_Add(
    "${NAME}_debug"
    GIT_REPOSITORY ${GIT_REPO}
    GIT_TAG ${GIT_TAG}
    PREFIX "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug"
    INSTALL_DIR "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug/install"
    CMAKE_ARGS
      -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ${ARGN}
  )
endfunction()

function(add_external_library NAME LIB_NAME)
  if(NOT NAME IN_LIST PRISMS_PF_VENDORED_PACKAGES)
    message(FATAL_ERROR "Invalid vendored package name")
    return()
  endif()

  add_library("imported_${NAME}" STATIC IMPORTED GLOBAL)
  add_library("imported_${NAME}_debug" STATIC IMPORTED GLOBAL)

  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/_deps/${NAME}/install/include")
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug/install/include")

  set_target_properties(
    "imported_${NAME}"
    PROPERTIES
      IMPORTED_LOCATION
        "${CMAKE_BINARY_DIR}/_deps/${NAME}/install/lib/${LIB_NAME}.a"
      INTERFACE_INCLUDE_DIRECTORIES
        "${CMAKE_BINARY_DIR}/_deps/${NAME}/install/include"
  )
  set_target_properties(
    "imported_${NAME}_debug"
    PROPERTIES
      IMPORTED_LOCATION
        "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug/install/lib/${LIB_NAME}.a"
      INTERFACE_INCLUDE_DIRECTORIES
        "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug/install/include"
  )

  add_dependencies("imported_${NAME}" "${NAME}")
  add_dependencies("imported_${NAME}_debug" "${NAME}_debug")
endfunction()

function(install_external_library NAME)
	if(NOT NAME IN_LIST PRISMS_PF_VENDORED_PACKAGES)
    message(FATAL_ERROR "Invalid vendored package name")
    return()
  endif()

  install(
    DIRECTORY
      "${CMAKE_BINARY_DIR}/_deps/${NAME}/install/"
    DESTINATION "${CMAKE_INSTALL_BINDIR}/../vendored/${NAME}"
  )
  install(
    DIRECTORY
      "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug/install/"
    DESTINATION "${CMAKE_INSTALL_BINDIR}/../vendored/${NAME}_debug"
  )
endfunction()

# Add the external projects
add_external_project(libassert https://github.com/jeremy-rifkin/libassert.git v2.2.1)

# Create the libraries
add_external_library(libassert libassert)

# Create install rules
install_external_library(libassert)

# Add libassert to the Release and Debug lists
prisms_pf_add_dependency_target(imported_libassert_debug DEBUG PUBLIC)
prisms_pf_add_dependency_target(imported_libassert RELEASE PUBLIC)
