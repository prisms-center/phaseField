#
# Macro that expands inst.in files for explicit template
# instantiations
#

macro(expand_template_instantiations _target _inst_in_files)
  set(_inst_targets "")
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

    # Get the directory for finding templates
    if(DEFINED PRISMS_PF_CORE_DIR)
      # We're in an application context, use the install directory
      set(_template_file "${PRISMS_PF_CORE_DIR}/lib/cmake/prisms_pf/templates")
    else()
      # We're in the main project context, use current source dir
      set(_template_file "${CMAKE_BINARY_DIR}/cmake/templates")
    endif()

    # Create a tmp inst file in case the command fails to
    # execute. This two level thing is necessary so that
    # we don't try and compile with incomplete instantiations
    set(_final_output "${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}")
    add_custom_command(
      OUTPUT
        "${_final_output}"
      DEPENDS
        ${_dependency}
        ${_template_file}
        ${CMAKE_CURRENT_SOURCE_DIR}/${_inst_in_file}
      COMMAND
        ${_command}
      ARGS
        ${_template_file} < 
        ${CMAKE_CURRENT_SOURCE_DIR}/${_inst_in_file} >
        "${_final_output}.tmp"
      COMMAND
        ${CMAKE_COMMAND}
      ARGS
        -E rename "${_final_output}.tmp"
        "${_final_output}"
    )

    # Append to a list of inst _inst_targets
    list(APPEND _inst_targets "${_final_output}")
  endforeach()

  # Create a target that depends on the generations
  # of template files to ensure that things are done
  # in the right order
  add_custom_target(${_target}_inst ALL DEPENDS ${_inst_targets})

  # Add this target again to the different build types if the override is off.
  if(NOT DEFINED BUILD_OVERRIDE)
    foreach(_build ${PRISMS_PF_BUILD_TYPES})
      # The build types are uppercase so we need to
      # convert them
      string(TOLOWER ${_build} _build_lowercase)

      # Add the target
      add_dependencies(${_target}_${_build_lowercase} ${_target}_inst)
    endforeach()
  else()
    # Add the target
    add_dependencies(${_target} ${_target}_inst)
  endif()
endmacro()
