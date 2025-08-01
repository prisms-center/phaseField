#!/bin/bash
# SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
# SPDX-License-Identifier: GNU Lesser General Public Version 2.1
#
# This script copies a file to all the applications in the tests/automatic_tests 
# directory and all the applications in the applications directory.
#
#
# Usage:
# ./contrib/utilities/copy_file_to_all_applications.sh <file_to_copy>
#

if test ! -d src -o ! -d include -o ! -d applications; then
	echo "This script must be run from the top-level directory of PRISMS-PF"
	exit 1
fi

# Check that a file was provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <file_to_copy>"
    exit 1
fi

file_to_copy=$1

# Check that the file exists
if [ ! -f "$file_to_copy" ]; then
    echo "File '$file_to_copy' does not exist."
    exit 1
fi

filename=$(basename "$file_to_copy")

# These are the required files for a valid application
required_files=("main.cc" "custom_pde.h")

copy_to_valid_apps_in_dir() {
    local parent_dir="$1"
    for app_dir in "$parent_dir"/*/; do
        [ -d "$app_dir" ] || continue

        valid=true
        for required_file in "${required_files[@]}"; do
            if [ ! -f "$app_dir/$required_file" ]; then
                valid=false
                break
            fi
        done

        if $valid; then
            cp "$file_to_copy" "$app_dir"
            echo "Copied to $app_dir"
        fi
    done
}

copy_to_valid_apps_in_dir "applications"
copy_to_valid_apps_in_dir "tests/automatic_tests"

echo "Done"
exit 0