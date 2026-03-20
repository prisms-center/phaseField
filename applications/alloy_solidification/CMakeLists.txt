##
#  CMake script for the PRISMS-PF applications
##

cmake_minimum_required(VERSION 3.18)

# Create a project for the application
project(
  myapp
  C
  CXX
)

# Find the PRISMS-PF package
find_package(
  prisms_pf
  QUIET
  HINTS
    ${PRISMS_PF_DIR}
    $ENV{PRISMS_PF_DIR}
)
if(PRISMS_PF_FOUND)
  # Include setup script
  include(${PRISMS_PF_CMAKE_DIR}/setup_application.cmake)

  # Set up the application
   # Enable compile commands export
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

  set(
    TARGET_SRC
    "${CMAKE_CURRENT_SOURCE_DIR}/main.cc"
  )

  # Create the executables
  foreach(_build ${APPLICATION_BUILD_TYPES})
    string(TOLOWER ${_build} _build_lowercase)
    string(TOUPPER ${_build} _build_uppercase)

    set(_target main-${_build_lowercase})

    add_executable(${_target} ${TARGET_SRC})

    # put the binary in the source folder
    set_target_properties(
      ${_target}
      PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY
          "${CMAKE_CURRENT_SOURCE_DIR}"
    )

    if(TARGET prisms_pf::prisms_pf_${_build_lowercase})
      target_link_libraries(${_target} PRIVATE prisms_pf::prisms_pf_${_build_lowercase})
      message(STATUS "Successfully linked ${_target} to prisms_pf_${_build_lowercase}")
    else()
      message(
        FATAL_ERROR
        "Target prisms_pf_${_build_lowercase} does not exist! "
        "Available build types: ${PRISMS_PF_BUILD_TYPES}"
      )
    endif()

    target_include_directories(
      ${_target}
      PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}
    )

    target_compile_options(
      ${_target}
      PRIVATE
        $<$<COMPILE_LANGUAGE:CXX>:
        ${PRISMS_PF_WARNING_FLAGS}
        ${PRISMS_PF_CXX_FLAGS}
        ${PRISMS_PF_CXX_FLAGS_${_build_uppercase}}>
    )
    target_link_options(
      ${_target}
      PRIVATE
        $<$<LINK_LANGUAGE:CXX>:
        ${PRISMS_PF_LINKER_FLAGS}
        ${PRISMS_PF_LINKER_FLAGS_${_build_uppercase}}>
    )

    set(BUILD_OVERRIDE "")
  endforeach()

else()
  message(FATAL_ERROR "PRISMS-PF package not found. Please make sure it is in path.")
endif()
