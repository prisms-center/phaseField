#
# Macro that expands inst.in files for explicit template
# instantiations
#

macro(expand_template_instantiations _target _inst_in_files)
  # Loop over the provided inst.in files
  foreach(_inst_in_file ${_inst_in_files})
    # Get the name for the output file. In other words,
    # just drop the .in bit
    string(REGEX REPLACE "\\.in$" "" _inst_file "${_inst_in_file}")

    if(CMAKE_CROSSCOMPILING)
      message(FATAL_ERROR "We don't support cross compiling")
    endif()

    # Set the command we need to execute to run 
    # expand_template_instantiations.cc 
    set(_command expand_template_instantiations_exe)
    set(_dependency expand_template_instantiations_exe)

    # Create a tmp inst file in case the command fails to 
    # execute. This two level thing is necessary so that
    # we don't try and compile with incomplete instantiations
    add_custom_command(
      OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}
      DEPENDS ${_dependency}
              ${CMAKE_BINARY_DIR}/cmake/templates
              ${CMAKE_CURRENT_SOURCE_DIR}/${_inst_in_file}
      COMMAND ${_command}
      ARGS ${CMAKE_BINARY_DIR}/cmake/templates
           < ${CMAKE_CURRENT_SOURCE_DIR}/${_inst_in_file}
           > ${CMAKE_CURRENT_SOURCE_DIR}/${_inst_file}.tmp
      COMMAND ${CMAKE_COMMAND}
      ARGS -E rename
           ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}.tmp
           ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}
    )

    # Append to a list of inst _inst_targets
    list(APPEND _inst_targets ${CMAKE_CURRENT_BINARY_DIR}/${_inst_file})
  endforeach()

  # Create a target that depends on the generations
  # of template files to ensure that things are done
  # in the right order
  add_custom_target(${_target}_inst ALL 
    DEPENDS ${_inst_targets}
  )

  # Add this target again to the different build types
  foreach(_build ${PRISMS_PF_BUILD_TYPES})
    # The build types are uppercase so we need to 
    # convert them
    string(TOLOWER ${_build} _build_lowercase)

    # Add the target
    add_dependencies(${_target}_${_build_lowercase} ${_target}_inst)
  endforeach()
endmacro()
