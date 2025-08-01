#
# Define a sublibrary of PRISMS-PF
#

function(define_library _library)
    foreach(_build ${PRISMS_PF_BUILD_TYPES})
        string(TOLOWER ${_build} _build_lowercase)
        string(TOUPPER ${_build} _build_uppercase)
        set(_target "${_library}_${_build_lowercase}")

        add_library(${_target} ${ARGN})

        set_target_properties(${_target} PROPERTIES LINKER_LANGUAGE "CXX")

        target_include_directories(
            ${_target}
            SYSTEM
            PUBLIC ${PRISMS_PF_INCLUDE_DIRS}
        )

        target_include_directories(
            ${_target}
            PRIVATE ${CMAKE_BINARY_DIR}/include ${CMAKE_BINARY_DIR}/src
        )

        if(${VTK_BUILT_SEPARATELY})
            include_directories(SYSTEM ${VTK_INCLUDE_DIR})
            target_link_libraries(${_target} ${VTK_NEW_LIBRARIES})
        endif()
        if(${PRISMS_PF_WITH_CALIPER})
            target_link_libraries(${_target} caliper)
        endif()

        deal_ii_setup_target(${_target} ${_build_uppercase})
    endforeach()
endfunction()
