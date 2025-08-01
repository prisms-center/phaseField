#
# Try and find the Caliper library
#

if(PRISMS_PF_AUTODETECTION)
    find_package(CALIPER QUIET HINTS ${CALIPER_DIR} $ENV{CALIPER_DIR})
    if(${CALIPER_FOUND})
        set(PRISMS_PF_WITH_CALIPER ON)
    endif()
endif()
if(PRISMS_PF_WITH_CALIPER)
    find_package(CALIPER QUIET HINTS ${CALIPER_DIR} $ENV{CALIPER_DIR})
    message(STATUS "CALIPER_DIR: ${CALIPER_DIR}")
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
