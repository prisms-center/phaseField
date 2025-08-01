#
# Get the git information
#

macro(prisms_pf_git_version)
    message(STATUS "Querying git information")

    find_package(Git QUIET)

    set(PRISMS_PF_GIT_BRANCH "")
    set(PRISMS_PF_GIT_REVISION "")
    set(PRISMS_PF_GIT_SHORTREV "")
    set(PRISMS_PF_GIT_TIMESTAMP "")

    set(tmp "")
    set(tmp1 "")
    set(tmp2 "")

    if(GIT_FOUND AND EXISTS ${CMAKE_SOURCE_DIR}/.git/HEAD)
        # Grab the nearest tag with some additional information
        execute_process(
            COMMAND git describe --abbrev=12 --dirty --always --tags
            WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
            OUTPUT_VARIABLE tmp
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )

        # Grab the branch
        execute_process(
            COMMAND git symbolic-ref HEAD
            WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
            OUTPUT_VARIABLE tmp1
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(NOT ${tmp1} STREQUAL "")
            string(
                REGEX REPLACE
                "refs/heads/"
                ""
                PRISMS_PF_GIT_BRANCH
                "${tmp1}"
            )
        endif()

        # Grab the latest commit hash and date
        set(date_formatting "")
        if(NOT ${GIT_VERSION_STRING} VERSION_LESS 2.2)
            set(date_formatting "--date=iso-strict")
        endif()
        execute_process(
            COMMAND
                git log -n 1
                --pretty=format:"revision=%H, shortrev=%h, date=%cd"
                ${date_formatting}
            WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
            OUTPUT_VARIABLE tmp2
            OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        if(NOT ${tmp2} STREQUAL "")
            string(
                REGEX REPLACE
                "^\"revision=(.+), shortrev=(.+), date=(.+)\"$"
                "\\1"
                PRISMS_PF_GIT_REVISION
                "${tmp2}"
            )
            string(
                REGEX REPLACE
                "^\"revision=(.+), shortrev=(.+), date=(.+)\"$"
                "\\2"
                PRISMS_PF_GIT_SHORTREV
                "${tmp2}"
            )
            string(
                REGEX REPLACE
                "^\"revision=(.+), shortrev=(.+), date=(.+)\"$"
                "\\3"
                PRISMS_PF_GIT_TIMESTAMP
                "${tmp2}"
            )
            string(
                REPLACE
                "T"
                " "
                PRISMS_PF_GIT_TIMESTAMP
                "${PRISMS_PF_GIT_TIMESTAMP}"
            )
        endif()
    else()
        message(WARNING "git repo not found")
    endif()

    set(PRISMS_PF_GIT_TAG "${tmp}" CACHE INTERNAL "")
    unset(tmp)
    unset(tmp1)
    unset(tmp2)
endmacro()
