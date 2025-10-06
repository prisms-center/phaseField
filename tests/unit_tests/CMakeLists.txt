##
# CMake for PRISMS-PF unit tests
##

message(STATUS "Setting up PRISMS-PF unit tests...")

# Add each test subdirectory
add_subdirectory(core)
add_subdirectory(field_input)
add_subdirectory(user_inputs)
add_subdirectory(utilities)

message(STATUS "PRISMS-PF unit tests setup complete.")
