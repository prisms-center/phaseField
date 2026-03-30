#
# Print variable name and value
#
function(print_variable VARIABLE)
  message(STATUS "${VARIABLE} = ${${VARIABLE}}")
endfunction()

#
# Print text in subsection format
#
function(print_subsection TEXT)
  message(STATUS "")
  message(STATUS "=========================================================")
  message(STATUS "${TEXT}")
  message(STATUS "=========================================================")
  message(STATUS "")
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
