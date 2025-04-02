#!/bin/bash
# SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
# SPDX-License-Identifier: GNU Lesser General Public Version 2.1
#
# This script runs cppcheck on the PRISMS-PF code base.
#
#
# Usage:
# ./contrib/utilities/cppcheck.sh CMAKE_ARGS
#   with:
#     CMAKE_ARGS are the optional arguments passed to cmake
#

if test ! -d src -o ! -d include -o ! -d applications; then
  echo "This script must be run from the top-level directory of PRISMS-PF"
  exit 1
fi

if ! [ -x "$(command -v cppcheck)" ]; then
  echo "make sure cppcheck is in the path"
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

# Create a file that contains all the headers to ensure that cppcheck runs on all headers
(
  cd include
  find . -name '*.h'
) | grep -v allheaders.h | sed 's|^./|#include <|' | sed 's|$|>|' >include/prismspf/allheaders.h

# Run cppcheck
cppcheck \
  --enable=all -q --language=c++ --std=c++17 --inline-suppr \
  --suppress=missingIncludeSystem --suppress=unknownMacro \
  --project=compile_commands.json --check-level=exhaustive \
  >"output.txt" 2>&1

# Separate unusedFunction messages into a different file. This is because most functions
# are used by the user. It doesn't hurt to periodically check.
grep 'unusedFunction' output.txt | sort | uniq >cppcheck_unused.log

# grep interesting errors and make sure we remove duplicates:
grep -E '(warning|error|style|performance|portability): ' output.txt | grep -v 'unusedFunction' | sort | uniq >cppcheck.log

# If we have errors, report them and set exit status to failure
if [ -s cppcheck.log ]; then
  cat cppcheck.log
  exit 4
fi

# Print unused functions
cat cppcheck_unused.log

echo "OK"
exit 0
