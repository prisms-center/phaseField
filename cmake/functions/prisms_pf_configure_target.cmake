#
# Add dependency target to lists
#
function(prisms_pf_add_dependency_target TARGET BUILD TYPE)
  string(TOUPPER "${BUILD}" BUILD)
  string(TOUPPER "${TYPE}" TYPE)

  set(LIST_NAME "PRISMS_PF_${TYPE}_PACKAGES_${BUILD}")

  list(APPEND ${LIST_NAME} ${TARGET})

  set(${LIST_NAME} ${${LIST_NAME}} PARENT_SCOPE)
endfunction()

#
# Configure targets
#
function(prisms_pf_configure_targets TARGETS)
  message(STATUS ${PRISMS_PF_PUBLIC_PACKAGES_RELEASE})
  message(STATUS ${PRISMS_PF_PUBLIC_PACKAGES_DEBUG})

  foreach(_target IN LISTS TARGETS)
    # Determine whether to use release or debug things
    set(_use_debug OFF)
    if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
      set(_use_debug ON)
    elseif(${CMAKE_BUILD_TYPE} STREQUAL "DebugRelease")
      if(${_target} MATCHES "_debug$")
        set(_use_debug ON)
      endif()
    endif()

    # Set the minimum C++ standard. This should be private so
    # we don't infect downstream targets
    # TODO: These shouldn't be public
    target_compile_definitions(${_target} PUBLIC cxx_std_20)
    set_target_properties(
      ${_target}
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
    if(_use_debug)
      target_compile_options(
        ${_target}
        PUBLIC
          ${PRISMS_PF_CXX_FLAGS_LIST}
          ${PRISMS_PF_CXX_FLAGS_DEBUG_LIST}
      )
    else()
      target_compile_options(
        ${_target}
        PUBLIC
          ${PRISMS_PF_CXX_FLAGS_LIST}
          ${PRISMS_PF_CXX_FLAGS_RELEASE_LIST}
      )
    endif()

    if(_use_debug)
      target_link_options(
        ${_target}
        PUBLIC
          ${DEAL_II_LINKER_FLAGS_LIST}
          ${DEAL_II_LINKER_FLAGS_DEBUG_LIST}
      )
    else()
      target_link_options(
        ${_target}
        PUBLIC
          ${DEAL_II_LINKER_FLAGS_LIST}
          ${DEAL_II_LINKER_FLAGS_RELEASE_LIST}
      )
    endif()

    # clang-tidy if defined
    if(CLANG_TIDY_TOOL)
      set_target_properties(
        ${_target}
        PROPERTIES
          CXX_CLANG_TIDY
            "${CLANG_TIDY_TOOL}"
      )
    endif()

    # Include the config file in the target
    target_include_directories(
      ${_target}
      PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )

    # Link dependencies
    if(_use_debug)
      target_link_libraries(
        ${_target}
        PUBLIC
          ${PRISMS_PF_PUBLIC_PACKAGES_DEBUG}
        PRIVATE
          ${PRISMS_PF_PRIVATE_PACKAGES_DEBUG}
        INTERFACE
          ${PRISMS_PF_INTERFACE_PACKAGES_DEBUG}
      )
    else()
      target_link_libraries(
        ${_target}
        PUBLIC
          ${PRISMS_PF_PUBLIC_PACKAGES_RELEASE}
        PRIVATE
          ${PRISMS_PF_PRIVATE_PACKAGES_RELEASE}
        INTERFACE
          ${PRISMS_PF_INTERFACE_PACKAGES_RELEASE}
      )
    endif()
  endforeach()
endfunction()
