#
# Autopilot function that will set up an application
#

macro(prisms_pf_autopilot PRISMS_PF_CORE_DIR)
  # Enable compile commands export
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

  include_directories(${CMAKE_BINARY_DIR})

  # Set location of files
  include_directories(${PRISMS_PF_CORE_DIR}/include)
  include_directories(${PRISMS_PF_CORE_DIR}/src)
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})

  # Find Caliper if enabled
  if(${PRISMS_PF_WITH_CALIPER})
    find_package(CALIPER REQUIRED)
    include_directories(${CALIPER_INCLUDE_DIR})
  endif()

  # Set the location of application source files
  # Check if an override is specified otherwise 
  # proceed with defaults
  if(NOT DEFINED TARGET_SRC_OVERRIDE)
    set(TARGET_SRC 
      "${CMAKE_CURRENT_SOURCE_DIR}/main.cc" 
      "${CMAKE_CURRENT_SOURCE_DIR}/equations.cc" 
      "${CMAKE_CURRENT_SOURCE_DIR}/ICs_and_BCs.cc"
    )
  else()
    set(TARGET_SRC ${TARGET_SRC_OVERRIDE})
  endif()

  # Set targets & link libraries for the build type
  if(${PRISMS_PF_BUILD_DEBUG} STREQUAL "ON")
    add_executable(main_debug ${TARGET_SRC})
    set_property(TARGET main_debug PROPERTY OUTPUT_NAME main-debug)
    deal_ii_setup_target(main_debug DEBUG)
    target_link_libraries(main_debug ${PRISMS_PF_CORE_DIR}/libprisms-pf-debug.a)

    if(${PRISMS_PF_WITH_CALIPER})
      target_link_libraries(main_debug caliper)
    endif()
  endif()

  if(${PRISMS_PF_BUILD_RELEASE} STREQUAL "ON")
    add_executable(main_release ${TARGET_SRC})
    set_property(TARGET main_release PROPERTY OUTPUT_NAME main)
    deal_ii_setup_target(main_release RELEASE)
    target_link_libraries(main_release ${PRISMS_PF_CORE_DIR}/libprisms-pf-release.a)

    if(${PRISMS_PF_WITH_CALIPER})
      target_link_libraries(main_release caliper)
    endif()
  endif()
endmacro() 
