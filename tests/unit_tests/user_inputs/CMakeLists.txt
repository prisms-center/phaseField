file(GLOB TEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")

add_executable(test_user_inputs ${TEST_SOURCES})

target_link_libraries(test_user_inputs PRIVATE prisms_pf_debug)

target_include_directories(test_user_inputs PRIVATE
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_BINARY_DIR}/include
    ${CMAKE_SOURCE_DIR}/contrib/catch
    ${CMAKE_CURRENT_SOURCE_DIR}
)

add_test(NAME user_inputs COMMAND test_user_inputs)