#!/bin/bash
# SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
# SPDX-License-Identifier: GNU Lesser General Public Version 2.1
#
# This script runs cppcheck on the PRISMS-PF code base.
#
#
# Usage:
# ./contrib/utilities/cppcheck.sh
#

if test ! -d src -o ! -d include -o ! -d applications; then
  echo "This script must be run from the top-level directory of PRISMS-PF"
  exit 1
fi

if ! [ -x "$(command -v cppcheck)" ]; then
  echo "make sure cppcheck is in the path"
  exit 1
fi

# Run cppcheck
# TODO: Check the applications
cppcheck \
  --enable=all --language=c++ --std=c++17 \
  --suppress=missingIncludeSystem --suppress=unknownMacro \
  include src tests/automatic_tests tests/performance_tests >"output.txt" 2>&1

# grep interesting errors and make sure we remove duplicates:
grep -E '(warning|error|style|performance|portability): ' output.txt | sort | uniq >cppcheck.log

# If we have errors, report them and set exit status to failure
if [ -s cppcheck.log ]; then
  cat cppcheck.log
  exit 2
fi

echo "OK"
exit 0
