file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")

add_executable(test_field_input ${TEST_SOURCES})

target_link_libraries(test_field_input PRIVATE prisms_pf_debug)

target_include_directories(
  test_field_input
  PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/include
    ${CMAKE_SOURCE_DIR}/contrib/catch
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(NAME field_input COMMAND test_field_input)
