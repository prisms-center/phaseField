#
# Add an external project that works with DebugRelease
#
function(prisms_pf_add_external_project NAME GIT_REPO GIT_TAG)
  if(NOT NAME IN_LIST PRISMS_PF_VENDORED_PACKAGES)
    message(FATAL_ERROR "Invalid vendored package name")
    return()
  endif()

  set(LIB_NAMES ${ARGN})
  set(BUILD_BYPRODUCTS_ARGS "")
  foreach(LIB_NAME IN LISTS LIB_NAMES)
    list(
      APPEND BUILD_BYPRODUCTS_ARGS
      "<INSTALL_DIR>/${CMAKE_INSTALL_LIBDIR}/${LIB_NAME}.a"
    )
  endforeach()

  if(CMAKE_BUILD_TYPE STREQUAL "DebugRelease")
    ExternalProject_Add(
      ${NAME}
      GIT_REPOSITORY ${GIT_REPO}
      GIT_TAG ${GIT_TAG}
      PREFIX "${CMAKE_BINARY_DIR}/_deps/${NAME}"
      INSTALL_DIR "${CMAKE_BINARY_DIR}/_deps/${NAME}"
      CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ${ARGN}
      BUILD_BYPRODUCTS
        ${BUILD_BYPRODUCTS_ARGS}
    )
    ExternalProject_Add(
      "${NAME}_debug"
      GIT_REPOSITORY ${GIT_REPO}
      GIT_TAG ${GIT_TAG}
      PREFIX "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug"
      INSTALL_DIR "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug"
      CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=Debug -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR> ${ARGN}
      BUILD_BYPRODUCTS
        ${BUILD_BYPRODUCTS_ARGS}
    )
  else()
    ExternalProject_Add(
      ${NAME}
      GIT_REPOSITORY ${GIT_REPO}
      GIT_TAG ${GIT_TAG}
      PREFIX "${CMAKE_BINARY_DIR}/_deps/${NAME}"
      INSTALL_DIR "${CMAKE_BINARY_DIR}/_deps/${NAME}"
      CMAKE_ARGS
        -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} -DCMAKE_INSTALL_PREFIX=<INSTALL_DIR>
        ${ARGN}
      BUILD_BYPRODUCTS
        ${BUILD_BYPRODUCTS_ARGS}
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
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/_deps/${NAME}/include")
    set_target_properties(
      "imported_${NAME}"
      PROPERTIES
        IMPORTED_LOCATION
          "${CMAKE_BINARY_DIR}/_deps/${NAME}/${CMAKE_INSTALL_LIBDIR}/${LIB_NAME}.a"
        INTERFACE_INCLUDE_DIRECTORIES
          "${CMAKE_BINARY_DIR}/_deps/${NAME}/include"
    )
    add_dependencies("imported_${NAME}" "${NAME}")

    add_library("imported_${NAME}_debug" STATIC IMPORTED GLOBAL)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug/include")
    set_target_properties(
      "imported_${NAME}_debug"
      PROPERTIES
        IMPORTED_LOCATION
          "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug/${CMAKE_INSTALL_LIBDIR}/${LIB_NAME}.a"
        INTERFACE_INCLUDE_DIRECTORIES
          "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug/include"
    )
    add_dependencies("imported_${NAME}_debug" "${NAME}_debug")
  else()
    add_library("imported_${NAME}" STATIC IMPORTED GLOBAL)
    file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/_deps/${NAME}/include")
    set_target_properties(
      "imported_${NAME}"
      PROPERTIES
        IMPORTED_LOCATION
          "${CMAKE_BINARY_DIR}/_deps/${NAME}/${CMAKE_INSTALL_LIBDIR}/${LIB_NAME}.a"
        INTERFACE_INCLUDE_DIRECTORIES
          "${CMAKE_BINARY_DIR}/_deps/${NAME}/include"
    )
    add_dependencies("imported_${NAME}" "${NAME}")
  endif()
endfunction()

#
# Add imported library for transitive dependencies that works with DebugRelease
#
# NOTE: This is meant for dependencies that imported by other projects. For example,
# libassert relies on cpptrace and a few other packages. To ensure that the links
# are intact for our unit tests we have to use this function.
#
function(prisms_pf_import_archive TARGET_NAME PREFIX_DIR ARCHIVE_NAME)
  if(CMAKE_BUILD_TYPE STREQUAL "DebugRelease")
    if(TARGET "${TARGET_NAME}_debug")
      return()
    endif()
    add_library("${TARGET_NAME}_debug" STATIC IMPORTED GLOBAL)
    set_target_properties(
      "${TARGET_NAME}_debug"
      PROPERTIES
        IMPORTED_LOCATION
          "${PREFIX_DIR}/${CMAKE_INSTALL_LIBDIR}/${ARCHIVE_NAME}.a"
    )
    if(IS_DIRECTORY "${PREFIX_DIR}_debug/include")
      set_target_properties(
        "${TARGET_NAME}_debug"
        PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES
            "${PREFIX_DIR}_debug/include"
      )
    endif()
  endif()

  if(TARGET ${TARGET_NAME})
    return()
  endif()
  add_library(${TARGET_NAME} STATIC IMPORTED GLOBAL)
  set_target_properties(
    ${TARGET_NAME}
    PROPERTIES
      IMPORTED_LOCATION
        "${PREFIX_DIR}/${CMAKE_INSTALL_LIBDIR}/${ARCHIVE_NAME}.a"
  )
  if(IS_DIRECTORY "${PREFIX_DIR}/include")
    set_target_properties(
      ${TARGET_NAME}
      PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES
          "${PREFIX_DIR}/include"
    )
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
        "${CMAKE_BINARY_DIR}/_deps/${NAME}/"
      DESTINATION "${CMAKE_INSTALL_BINDIR}/../vendored/${NAME}"
      PATTERN "src"
        EXCLUDE
      PATTERN "tmp"
        EXCLUDE
    )
    install(
      DIRECTORY
        "${CMAKE_BINARY_DIR}/_deps/${NAME}_debug/"
      DESTINATION "${CMAKE_INSTALL_BINDIR}/../vendored/${NAME}_debug"
      PATTERN "src"
        EXCLUDE
      PATTERN "tmp"
        EXCLUDE
    )
  else()
    install(
      DIRECTORY
        "${CMAKE_BINARY_DIR}/_deps/${NAME}/"
      DESTINATION "${CMAKE_INSTALL_BINDIR}/../vendored/${NAME}"
      PATTERN "src"
        EXCLUDE
      PATTERN "tmp"
        EXCLUDE
    )
  endif()
endfunction()
