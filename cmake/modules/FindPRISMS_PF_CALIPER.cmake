#
# Try and find the Caliper library
#

set(CALIPER_DIR "" CACHE PATH "An optional hint to a Caliper directory")
set_if_empty(CALIPER_DIR "$ENV{CALIPER_DIR}")

# As an optional feature, we must either let the user specify to use Caliper
# explicitly, or we can try and auto-detect it
if(PRISMS_PF_AUTODETECTION)
  message(STATUS "Attempting to auto-detect Caliper installation...")
  message(STATUS "  CALIPER_DIR: ${CALIPER_DIR}")
  find_package(Caliper QUIET HINTS ${CALIPER_DIR})
  if(${CALIPER_FOUND})
    set(PRISMS_PF_WITH_CALIPER ON)
    message(STATUS "  Caliper found at ${CALIPER_DIR}")
  else()
    set(PRISMS_PF_WITH_CALIPER OFF)
    message(STATUS "  Caliper not found")
  endif()
endif()
if(PRISMS_PF_WITH_CALIPER)
  message(STATUS "Looking for Caliper installation...")
  message(STATUS "  CALIPER_DIR: ${CALIPER_DIR}")
  find_package(Caliper QUIET HINTS ${CALIPER_DIR})
  if(${CALIPER_FOUND})
    include_directories(${CALIPER_INCLUDE_DIR})
    message(STATUS "  Caliper found at ${CALIPER_DIR}")
  else()
    message(
      FATAL_ERROR
      "Caliper not found. Disable PRISMS_PF_WITH_CALIPER or specify a hint to your installation directory with CALIPER_DIR"
    )
  endif()
endif()
