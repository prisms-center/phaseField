#
# Wrapper for the target_sources function.
#
# This takes a list of targets, header, source, and inst
# files and add them to the given target, making sure to
# add code generator dependencies for the inst files.
#

function(prisms_pf_target_sources targets headers sources inst_bases)
  # Derive the .inst.in list from the base names
  set(_inst)
  foreach(_base ${inst_bases})
    list(APPEND _inst "${_base}.inst.in")
  endforeach()

  # Use the first target to drive .inst generation.
  # We don't want this to happen twice.
  list(GET targets 0 _primary_target)
  expand_template_instantiations(${_primary_target} "${_inst}")

  # Add sources and includes to each target
  foreach(_target IN LISTS targets)
    # Add sources to target
    target_sources(
      ${_target}
      PRIVATE
        ${sources}
      PUBLIC
        FILE_SET HEADERS
          BASE_DIRS ${PROJECT_SOURCE_DIR}/include
          FILES ${headers}
    )

    # Includes
    target_include_directories(${_target} PRIVATE ${CMAKE_BINARY_DIR}/src)

    # Let CMake know that .cc files can't compile until the .inst
    # files exist
    foreach(_base ${_inst_bases})
      set_source_files_properties(
        ${_base}.cc
        PROPERTIES
          OBJECT_DEPENDS
            "${CMAKE_CURRENT_BINARY_DIR}/${_base}.inst"
      )
    endforeach()
  endforeach()
endfunction()
