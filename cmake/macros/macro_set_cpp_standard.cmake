#
# Set the C++ standard
#

macro(set_cpp_standard _number)
  set(CMAKE_CXX_STANDARD ${_number})
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endmacro()
