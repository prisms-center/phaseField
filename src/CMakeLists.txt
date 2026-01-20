# This file collects the source files for PRISMS-PF
include_directories(${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/src)
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${CMAKE_BINARY_DIR}/include)

foreach(TARGET IN ITEMS core solvers user_inputs utilities)
  add_subdirectory(${TARGET})
endforeach()

# Handles the install path
foreach(_target ${PRISMS_PF_TARGETS})
  target_include_directories(${_target} PUBLIC $<INSTALL_INTERFACE:include>)
endforeach()

# Add the dummy.cc
set(_src dummy.cc)
set(_header)
collect_header_and_src_files("${_src}" "${_header}")
