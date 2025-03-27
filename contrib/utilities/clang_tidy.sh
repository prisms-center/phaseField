#!/bin/bash
# SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
# SPDX-License-Identifier: GNU Lesser General Public Version 2.1
#
# This script runs clang-tidy on the PRISMS-PF code base.
#
#
# Usage:
# ./contrib/utilities/clang_tidy.sh CMAKE_ARGS
#   with:
#     CMAKE_ARGS are the optional arguments passed to cmake
#

if test ! -d src -o ! -d include -o ! -d applications; then
    echo "This script must be run from the top-level directory of PRISMS-PF"
    exit 1
fi

# TODO: Don't think I need to use clang compiler if mpi compiler is being used by default
if [ -x "$(command -v run-clang-tidy)" ] && [ -x "$(command -v clang++)" ]; then
    CLANG_TIDY=run-clang-tidy
    CXX_COMPILER=clang++
elif [ -x "$(command -v run-clang-tidy-18)" ] && [ -x "$(command -v clang++-18)" ]; then
    CLANG_TIDY=run-clang-tidy-18
    CXX_COMPILER=clang++-18
else
    echo "Neither run-clang-tidy and clang++ nor run-clang-tidy-18 and clang++-18 are in path."
    exit 1
fi

# Construct the cmake arguments
ARGS=("-D" "CMAKE_EXPORT_COMPILE_COMMANDS=ON" "-D" "CMAKE_BUILD_TYPE=Debug" "$@")

# Compile
if [ -f CMakeCache.txt ]; then
    rm -f CMakeCache.txt
fi
cmake "${ARGS[@]}" . || (
    echo "cmake failed!"
    false
) || exit 2
cmake --build . || exit 3

# Create a file that contains all the headers to ensure that clang-tidy runs on all headers
(
    cd include
    find . -name '*.h'
) | grep -v allheaders.h | sed 's|^./|#include <|' | sed 's|$|>|' >include/allheaders.h

# TODO: cat the clang-tidy rules here

# Run clang-tidy
$CLANG_TIDY -p . -quiet -header-filter "include/*" -extra-arg='-DCLANG_TIDY' 2>error.txt >output.txt

# grep interesting errors and make sure we remove duplicates:
grep -E '(warning|error): ' output.txt | sort | uniq >clang-tidy.log

# If we have errors, report them and set exit status to failure
if [ -s clang-tidy.log ]; then
    cat clang-tidy.log
    exit 4
fi

echo "OK"
exit 0
