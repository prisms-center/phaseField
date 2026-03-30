#
# Configure a target given a build type.
#

function(prisms_pf_configure_target target_name build_type)
  # Set the minimum C++ standard. This should be private so
  # we don't infect downstream targets
  # TODO: These shouldn't be public
  target_compile_definitions(${target_name} PUBLIC cxx_std_20)
  set_target_properties(
    ${target_name}
    PROPERTIES
      CXX_EXTENSIONS
        OFF
      CXX_STANDARD
        20
      CXX_STANDARD_REQUIRED
        ON
  )

  # Add the compile flags, which we inherit from deal.II
  # TODO: We really shouldn't do this and have suggested
  # flags that get inherited from deal.II
  # TODO: These shouldn't be public
  target_compile_options(
    ${target_name}
    PUBLIC
      ${DEAL_II_CXX_FLAGS_LIST}
      $<$<CONFIG:Debug>:${DEAL_II_CXX_FLAGS_DEBUG_LIST}>
      $<$<CONFIG:Release>:${DEAL_II_CXX_FLAGS_RELEASE_LIST}>
  )

  target_link_options(
    ${target_name}
    PUBLIC
      SHELL:${DEAL_II_LINKER_FLAGS}
      $<$<CONFIG:Debug>:SHELL:${DEAL_II_LINKER_FLAGS_DEBUG}>
      $<$<CONFIG:Release>:SHELL:${DEAL_II_LINKER_FLAGS_RELEASE}>
  )

  # Include the config file in the target
  target_include_directories(
    ${target_name}
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
      $<INSTALL_INTERFACE:include>
  )

  # Link the matching deal.II target

  # Again we run into some more of the uniqueness of deal.II with the way they set
  # up build types. Their default build type is `DebugRelease`, which builds both
  # Debug and Release targets. As such, our project much choose what deal.II target
  # to link to when the user compiles with Debug and Release.
  #
  # When deal.II is built with only Debug or only Release, we link against for all
  # PRISMS-PF builds. For example, if deal.II is built in Release, PRISMS-PF in
  # Debug will only contain assertions and debugging for our library. In most cases,
  # this is sufficient to catch most issues.
  #
  # When deal.II is build with DebugRelease, we link to the target that corresponds
  # to the PRISMS-PF build type. This is recommended for developers.
  if(DEAL_II_BUILD_TYPE STREQUAL "DebugRelease")
    if(build_type STREQUAL "Debug")
      target_link_libraries(${target_name} PUBLIC dealii::dealii_debug)
    elseif(build_type STREQUAL "Release")
      target_link_libraries(${target_name} PUBLIC dealii::dealii_release)
    else()
      message(FATAL_ERROR "Unknown build_type = ${build_type}")
    endif()
  elseif(DEAL_II_BUILD_TYPE STREQUAL "Debug")
    target_link_libraries(${target_name} PUBLIC dealii::dealii_debug)
  elseif(DEAL_II_BUILD_TYPE STREQUAL "Release")
    target_link_libraries(${target_name} PUBLIC dealii::dealii_release)
  else()
    message(FATAL_ERROR "Unknown DEAL_II_BUILD_TYPE = ${DEAL_II_BUILD_TYPE}")
  endif()

  # Link other dependencies
  # TODO: We don't need to set everything as public, and we shouldn't.
  # Ideally, the only public bits need to be MPI and deal.II
  target_link_libraries(
    ${target_name}
    PUBLIC
      MPI::MPI_CXX
    PRIVATE
      $<BUILD_LOCAL_INTERFACE:libassert::assert>
  )
endfunction()
