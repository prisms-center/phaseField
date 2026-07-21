#
# Add an external project that works with DebugRelease
#
function(prisms_pf_add_external_project NAME GIT_REPO GIT_TAG)
  if(NOT NAME IN_LIST PRISMS_PF_VENDORED_PACKAGES)
    message(FATAL_ERROR "Invalid vendored package name")
    return()
  endif()

  if(CMAKE_BUILD_TYPE STREQUAL "DebugRelease")
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
  else()
    ExternalProject_Add(
      ${NAME}
      GIT_REPOSITORY ${GIT_REPO}
      GIT_TAG ${GIT_TAG}
      PREFIX "${CMAKE_BINARY_DIR}/_deps/${NAME}"
      INSTALL_DIR "${CMAKE_BINARY_DIR}/_deps/${NAME}/install"
      CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        ${ARGN}
    )
  endif()
endfunction()

#
# Add an external project library that works with DebugRelease
#
function(prisms_pf_add_external_library NAME LIB_NAME)
  if(NOT NAME IN_LIST PRISMS_PF_VENDORED_PACKAGES)
    message(FATAL_ERROR "Invalid vendored package name")
    return()
  endif()

  if(CMAKE_BUILD_TYPE STREQUAL "DebugRelease")
    add_library("imported_${NAME}" STATIC IMPORTED GLOBAL)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/_deps/${NAME}/install/include")
    set_target_properties(
      "imported_${NAME}"
      PROPERTIES
        IMPORTED_LOCATION
          "${CMAKE_BINARY_DIR}/_deps/${NAME}/install/lib/${LIB_NAME}.a"
        INTERFACE_INCLUDE_DIRECTORIES
          "${CMAKE_BINARY_DIR}/_deps/${NAME}/install/include"
    )
    add_dependencies("imported_${NAME}" "${NAME}")

    add_library("imported_${NAME}_debug" STATIC IMPORTED GLOBAL)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug/install/include")
    set_target_properties(
      "imported_${NAME}_debug"
      PROPERTIES
        IMPORTED_LOCATION
          "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug/install/lib/${LIB_NAME}.a"
        INTERFACE_INCLUDE_DIRECTORIES
          "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug/install/include"
    )
    add_dependencies("imported_${NAME}_debug" "${NAME}_debug")
  else()
    add_library("imported_${NAME}" STATIC IMPORTED GLOBAL)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/_deps/${NAME}/install/include")
    set_target_properties(
      "imported_${NAME}"
      PROPERTIES
        IMPORTED_LOCATION
          "${CMAKE_BINARY_DIR}/_deps/${NAME}/install/lib/${LIB_NAME}.a"
        INTERFACE_INCLUDE_DIRECTORIES
          "${CMAKE_BINARY_DIR}/_deps/${NAME}/install/include"
    )
    add_dependencies("imported_${NAME}" "${NAME}")
  endif()
endfunction()

#
# Add an external project install targets that works with DebugRelease
#
function(prisms_pf_install_external_library NAME)
  if(NOT NAME IN_LIST PRISMS_PF_VENDORED_PACKAGES)
    message(FATAL_ERROR "Invalid vendored package name")
    return()
  endif()

  if(CMAKE_BUILD_TYPE STREQUAL "DebugRelease")
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
  else()
    install(
      DIRECTORY
        "${CMAKE_BINARY_DIR}/_deps/${NAME}/install/"
      DESTINATION "${CMAKE_INSTALL_BINDIR}/../vendored/${NAME}"
    )
  endif()
endfunction()
