##
#  CMake script for the PRISMS-PF doxygen documentation.
##

cmake_minimum_required(VERSION 3.17.0)

project(PRISMS-PF)

# Check that Doxygen is installed
find_package(Doxygen REQUIRED)

if(DOXYGEN_FOUND)
    # Path of Doxyfile
    set(DOXYGEN_INPUT ${CMAKE_SOURCE_DIR}/Doxyfile)

    # Path of output file
    set(DOXYGEN_OUTPUT ${CMAKE_SOURCE_DIR}/Doxyfile.out)

    # Copy input to output
    configure_file(${DOXYGEN_INPUT} ${DOXYGEN_OUTPUT} @ONLY)

    # Add target so we can just call make with no target
    add_custom_target(
        doc_doxygen
        ALL
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUTPUT}
        COMMENT "Generating Doxygen documentation"
        VERBATIM
    )

    message(
        STATUS
        "Doxygen documentation will be generated using: ${DOXYGEN_INPUT}"
    )
endif(DOXYGEN_FOUND)
