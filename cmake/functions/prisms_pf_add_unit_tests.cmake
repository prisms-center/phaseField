#
# Add unit tests belonging to some module
#

function(add_unit_tests MODULE)
  set(_unit_tests ${ARGN})

  foreach(_test ${_unit_tests})
    string(REGEX REPLACE "\\.cc$" "" _test_name "${_test}")

    add_executable(${_test_name} ${_test})

    target_compile_features(${_test_name} PRIVATE cxx_std_20)
    set_target_properties(
      ${_test_name}
      PROPERTIES
        CXX_EXTENSIONS
          OFF
    )

    target_link_libraries(
      ${_test_name}
      PRIVATE
        prisms_pf::prisms_pf
        Catch2::Catch2WithMain
    )

    catch_discover_tests(${_test_name}
      TEST_PREFIX "${MODULE}::"
    )
  endforeach()
endfunction()
