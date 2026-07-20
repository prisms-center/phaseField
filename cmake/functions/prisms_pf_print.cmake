#
# Format variable into string and value
#
function(format_variable OUTPUT VARIABLE)
  set(${OUTPUT} "${VARIABLE} = ${${VARIABLE}}\n" PARENT_SCOPE)
endfunction()
#
# Print variable name and value
#
function(print_variable VARIABLE)
  format_variable(FORMATTED "${VARIABLE}")
  message(STATUS "${FORMATTED}")
endfunction()

#
# Format text as a subsection
#
function(format_subsection OUTPUT TEXT)
  set(
    ${OUTPUT}
    "\n=========================================================\n${TEXT}\n=========================================================\n"
    PARENT_SCOPE
  )
endfunction()

#
# Print text in subsection format
#
function(print_subsection TEXT)
  format_subsection(FORMATTED "${TEXT}")
  message(STATUS "${FORMATTED}")
endfunction()

#
# Print list of variable names and value with alignment
#
# NOTE: This can take blank strings
#
function(print_variables)
  # Find the longest variable name so we can align them
  set(MAX_LEN 0)
  math(EXPR LAST "${ARGC} - 1")
  foreach(I RANGE ${LAST})
    string(LENGTH "${ARGV${I}}" LEN)
    if(LEN GREATER MAX_LEN)
      set(MAX_LEN ${LEN})
    endif()
  endforeach()

  # Print each variable padded to align the values
  foreach(I RANGE ${LAST})
    set(VARIABLE "${ARGV${I}}")
    if(VARIABLE STREQUAL "")
      message(STATUS "")
    else()
      string(LENGTH "${VARIABLE}" LEN)
      math(EXPR PAD_LEN "${MAX_LEN} - ${LEN}")
      string(
        REPEAT " "
        ${PAD_LEN}
        PADDING
      )
      message(STATUS "${VARIABLE}:${PADDING} ${${VARIABLE}}")
    endif()
  endforeach()
endfunction()

#
# Print text to summary.log
#
function(print_to_summary)
  file(APPEND "${CMAKE_BINARY_DIR}/summary.log" "${ARGN}")
endfunction()

#
# Print text to detailed.log
#
function(print_to_detailed)
  file(APPEND "${CMAKE_BINARY_DIR}/detailed.log" "${ARGN}")
endfunction()

#
# Print text to both summary.log and detailed.log
#
function(print_to_both)
  print_to_summary("${ARGN}")
  print_to_detailed("${ARGN}")
endfunction()
