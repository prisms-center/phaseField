#
# Write the cmake config to a file
#

# Set the file pathes and remove any log files that may already exist
set(detailed_log "${CMAKE_SOURCE_DIR}/detailed.log")
set(summary_log "${CMAKE_SOURCE_DIR}/summary.log")
file(REMOVE ${detailed_log} ${summary_log})

function(write_to_both)
    file(APPEND "${detailed_log}" "${ARGN}")
    file(APPEND "${summary_log}" "${ARGN}")
endfunction()

function(write_to_summary)
    file(APPEND "${summary_log}" "${ARGN}")
endfunction()

function(write_to_detailed)
    file(APPEND "${detailed_log}" "${ARGN}")
endfunction()

write_to_both(
"#
#  PRISMS-PF configuration:
#    Version:               ${PRISMS_PF_VERSION}
#    Git tag:               ${PRISMS_PF_GIT_TAG}
#    Git branch:            ${PRISMS_PF_GIT_BRANCH}
#    Git short revision:    ${PRISMS_PF_GIT_SHORTREV}
"
)
write_to_detailed("#    Git revision:          ${PRISMS_PF_GIT_REVISION}\n")
write_to_both("#    Git timestamp:         ${PRISMS_PF_GIT_TIMESTAMP}
#
#    CMAKE_BUILD_TYPE:      ${CMAKE_BUILD_TYPE}
#    CMAKE_INSTALL_PREFIX:  ${CMAKE_INSTALL_PREFIX}
#    CMAKE_SOURCE_DIR:      ${CMAKE_SOURCE_DIR}
#    CMAKE_BINARY_DIR:      ${CMAKE_BINARY_DIR}
#    CMAKE_CXX_COMPILER:    ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} on platform ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}
#                           ${CMAKE_CXX_COMPILER}
"
)
write_to_detailed("#    CMAKE_GENERATOR:       ${CMAKE_GENERATOR}\n"
)
write_to_both("#
#    C++ standard:          ${CMAKE_CXX_STANDARD}
"
)
write_to_detailed("#
#    CXX_FLAGS:             ${DEAL_II_CXX_FLAGS}
#    CXX_FLAGS_DEBUG:       ${DEAL_II_CXX_FLAGS_DEBUG}
#    CXX_FLAGS_RELEASE:     ${DEAL_II_CXX_FLAGS_RELEASE}
#    ADDITIONAL_CXX_FLAGS:  ${PRISMS_PF_ADDITIONAL_CXX_FLAGS}\n"
)
write_to_both("#
#    DEAL_II_DIR:           ${deal.II_DIR}
#    DEAL_II_VERSION:       ${DEAL_II_PACKAGE_VERSION}
#
#    PRISMS_PF_AUTODETECTION: ${PRISMS_PF_AUTODETECTION}
#
#    PRISMS_PF_WITH_ZLIB:      ${PRISMS_PF_WITH_ZLIB}
#    PRISMS_PF_WITH_HDF5:      ${PRISMS_PF_WITH_HDF5}
#    PRISMS_PF_WITH_SUNDIALS:  ${PRISMS_PF_WITH_SUNDIALS}
#    PRISMS_PF_WITH_CALIPER:   ${PRISMS_PF_WITH_CALIPER}
#    PRISMS_PF_WITH_CUDA:      ${PRISMS_PF_WITH_CUDA}
#
#    64BIT_INDICES:            ${64BIT_INDICES}
#    ADDITIONAL_OPTIMIZATIONS: ${ADDITIONAL_OPTIMIZATIONS}
#    ADDITIONAL_DEGREES:       ${ADDITIONAL_DEGREES}
#    UNWRAP_COMPILER:          ${UNWRAP_COMPILER}
#"
)
