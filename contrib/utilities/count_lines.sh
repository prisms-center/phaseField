#!/bin/bash
# SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
# SPDX-License-Identifier: GNU Lesser General Public Version 2.1
#
# This script counts the lines in the PRISMS-PF code base.
#
#
# Usage:
# ./count_lines.sh PATH_TO_PRISMS_PF
#   with:
#     PATH_TO_PRISMS_PF as the path to the top level of the PRISMS-PF repository
#
# The output of the script is a list of dates and the number of lines in the source
# code and tests.
#

if [ -z "$1" ]; then
    echo "Error: No path provided."
    echo "Usage: ./count_lines.sh PATH_TO_PRISMS_PF"
    exit 1
fi
cd "$1" || {
    echo "Failed to cd into '$1'"
    exit 1
}
if test ! -d src -o ! -d include -o ! -d applications; then
    echo "This script must be run from the top-level directory of PRISMS-PF"
    exit 1
fi

current_branch=$(git rev-parse --abbrev-ref HEAD)
git checkout -q devel

commits=$(git log |
    grep -E '^commit ' |
    perl -p -e 's/^commit //g;' |
    perl -e '$i=0; while (<>) { ++$i; if ($i % 200 == 0) { print; } }')

for commit in $commits; do
    git checkout -q "$commit"

    date=$(git log --date=short --format="%ad" -n 1)

    files_source=$(find . -name '*.h' -or -name '*.cc' -type f |
        grep -E -i -v '(tests)/')
    lines_source=$(cat $files_source | wc -l)

    files_tests=$(find . -name '*.h' -or -name '*.cc' -type f |
        grep tests/)
    lines_tests=$(cat $files_tests | wc -l)

    echo "$date" "$lines_source" "$lines_tests"
done

git checkout -q ${current_branch}
