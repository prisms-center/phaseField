# This file collects the source files for PRISMS-PF
foreach(TARGET IN ITEMS core solvers user_inputs utilities)
    add_subdirectory(${TARGET})
endforeach()

# Add the dummy.cc
set(_src dummy.cc)
set(_header)
collect_header_and_src_files("${_src}" "${header}")
