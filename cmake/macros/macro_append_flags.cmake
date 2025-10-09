#
# Split a string of flags into a list. For example, this would take
# something like `-Wall -Wpedantic` and turn it into `-Wall; -Wpedantic`
#

macro(append_flags SOURCE_FLAGS DEST_FLAGS)
    if(${SOURCE_FLAGS})
        separate_arguments(_temp_flags NATIVE_COMMAND "${${SOURCE_FLAGS}}")
        foreach(flag IN LISTS _temp_flags)
            list(APPEND ${DEST_FLAGS} "${flag}")
        endforeach()
    endif()
endmacro()