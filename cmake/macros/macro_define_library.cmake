#
# Define a sublibrary of PRISMS-PF
#

function(define_library _library)
    foreach(_build ${PRISMS_PF_BUILD_TYPES})
        # Create a target for each build type. Note that the build types are
        # mixed case (e.g., Debug) so we need to convert them to lowercase and
        # uppercase for various uses.
        string(TOLOWER ${_build} _build_lowercase)
        string(TOUPPER ${_build} _build_uppercase)
        set(_target "${_library}_${_build_lowercase}")

        # Add the library and set some properties
        add_library(${_target} ${ARGN})
        set_target_properties(${_target} PROPERTIES LINKER_LANGUAGE CXX)

        # Set the include directories
        target_include_directories(
            ${_target}
            SYSTEM
            PUBLIC ${PRISMS_PF_INCLUDE_DIRS}
            PRIVATE ${CMAKE_BINARY_DIR}/include ${CMAKE_BINARY_DIR}/src
        )

        target_compile_options(
            ${_target}
            PRIVATE
                $<$<COMPILE_LANGUAGE:CXX>:${PRISMS_PF_WARNING_FLAGS} ${PRISMS_PF_CXX_FLAGS} ${PRISMS_PF_CXX_FLAGS_${_build_uppercase}}>
        )
        target_link_options(
            ${_target}
            PRIVATE
                $<$<LINK_LANGUAGE:CXX>:${PRISMS_PF_LINKER_FLAGS} ${PRISMS_PF_LINKER_FLAGS_${_build_uppercase}}>
        )

        # Add other dependencies, making sure they are public so that they
        # propagate to targets that link against this library

        # VTK
        if(${VTK_BUILT_SEPARATELY})
            target_include_directories(
                ${_target}
                SYSTEM
                PUBLIC ${VTK_INCLUDE_DIR}
            )
            target_link_libraries(${_target} PUBLIC ${VTK_NEW_LIBRARIES})
        endif()

        # caliper
        if(${PRISMS_PF_WITH_CALIPER})
            target_link_libraries(${_target} PUBLIC caliper)
        endif()

        # deal.II
        target_link_libraries(
            ${_target}
            PUBLIC ${DEAL_II_TARGET_${_build_uppercase}}
        )
    endforeach()
endfunction()
