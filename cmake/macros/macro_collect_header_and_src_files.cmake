#
# Collect the header and source files to add to the PRISMS-PF variables
#

macro(collect_header_and_src_files _sources _headers)
  set(_processed_sources "")
  set(_processed_headers "")

  foreach(_file ${_sources})
    if(NOT IS_ABSOLUTE "${_file}")
      list(APPEND _processed_sources "${CMAKE_CURRENT_SOURCE_DIR}/${_file}")
    else()
      list(APPEND _processed_sources "${_file}")
    endif()
  endforeach()
  foreach(_file ${_headers})
    if(NOT IS_ABSOLUTE "${_file}")
      list(APPEND _processed_headers "${CMAKE_CURRENT_SOURCE_DIR}/${_file}")
    else()
      list(APPEND _processed_headers "${_file}")
    endif()
  endforeach()

  set_property(
    GLOBAL
    APPEND
    PROPERTY
      PRISMS_PF_SOURCE_FILES
        ${_processed_sources}
  )
  set_property(
    GLOBAL
    APPEND
    PROPERTY
      PRISMS_PF_HEADER_FILES
        ${_processed_headers}
  )
endmacro()
