#
# Expands .inst.in files into .inst files using the
# expand_template_instantiations.cc executable. This
# allows us to automatically generatre explicit
# template instantiations in our source files.
#

function(expand_template_instantiations TARGET INST_IN_FILES)
  # Empty variable for the generated inst files that we'll add
  # to a target
  set(_inst_outputs "")

  # The cmake script that handles generator execution
  set(_script "${CMAKE_SOURCE_DIR}/cmake/scripts/expand_inst.cmake")

  # The template lookup file that's in the build tree
  set(_templates "${CMAKE_BINARY_DIR}/cmake/templates")
  foreach(_inst_in_file ${INST_IN_FILES})
    # Drop the .in suffix for the output filename
    string(REGEX REPLACE "\\.in$" "" _inst_file "${_inst_in_file}")

    # Input files from source tree and output in the build tree
    set(_input "${CMAKE_CURRENT_SOURCE_DIR}/${_inst_in_file}")
    set(_output "${CMAKE_CURRENT_BINARY_DIR}/${_inst_file}")

    # Make the custom target name unique per call site using the current directory
    string(MD5 _dir_hash "${CMAKE_CURRENT_SOURCE_DIR}")
    set(_inst_target "${TARGET}_inst_${_dir_hash}")

    add_custom_command(
      OUTPUT
        "${_output}"
      DEPENDS
        expand_template_instantiations_exe
        "${_templates}"
        "${_input}"
      COMMAND
        "${CMAKE_COMMAND}" -D EXE=$<TARGET_FILE:expand_template_instantiations_exe> -D
        TEMPLATES=${_templates} -D INPUT=${_input} -D OUTPUT=${_output} -D
        OUTPUT_TMP=${_output}.tmp -P "${_script}"
      COMMENT "Expanding instantiations: ${_inst_file}"
      VERBATIM
    )

    # Append to the list of inst files
    list(APPEND _inst_outputs "${_output}")
  endforeach()

  # Named target so the main library can depend on it cleanly
  add_custom_target(${_inst_target} ALL DEPENDS ${_inst_outputs})
  add_dependencies(${TARGET} ${_inst_target})

  # Return outputs to parent scope
  set(INST_OUTPUTS ${_inst_outputs} PARENT_SCOPE)
endfunction()
