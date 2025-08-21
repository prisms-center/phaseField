##
#  CMake script for the PRISMS-PF doxygen documentation.
##

# cmake config
message(STATUS "Using CMake ${CMAKE_VERSION}")
message(STATUS "")
cmake_minimum_required(VERSION 3.13.4)
cmake_policy(VERSION 3.13.4)

# Grab the version of PRISMS-PF
file(STRINGS "${CMAKE_SOURCE_DIR}/../../VERSION" PRISMS_PF_VERSION LIMIT_COUNT 1)

message(STATUS "")
message(STATUS "=========================================================")
message(STATUS "Configuring PRISMS-PF v${PRISMS_PF_VERSION} Documentation")
message(STATUS "=========================================================")
message(STATUS "")

project(PRISMS-PF)

# Check that Doxygen is installed
find_package(Doxygen REQUIRED)

if(DOXYGEN_FOUND)
    # Path of Doxyfile
    set(DOXYGEN_INPUT ${CMAKE_SOURCE_DIR}/Doxyfile)

    # Path of output file
    set(DOXYGEN_OUTPUT ${CMAKE_SOURCE_DIR}/Doxyfile.out)

    # The directory where the documentation will be generated
    set(DOXYGEN_HTML_OUTPUT "${PRISMS_PF_VERSION}")
    set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/docs")

    # Copy input to output
    configure_file(${DOXYGEN_INPUT} ${DOXYGEN_OUTPUT} @ONLY)

    # Add target to generate documentation
    doxygen_add_docs(
        doc
        ${CMAKE_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/include
        COMMENT "Generating Doxygen documentation"
        CONFIG_FILE ${DOXYGEN_OUTPUT}
        ALL
    )

    message(
        STATUS
        "Doxygen documentation will be generated using: ${DOXYGEN_INPUT}"
        "and output to: ${DOXYGEN_OUTPUT_DIRECTORY}/${DOXYGEN_HTML_OUTPUT}"
    )
endif(DOXYGEN_FOUND)
