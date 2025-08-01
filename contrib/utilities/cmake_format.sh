#!/bin/bash
# SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
# SPDX-License-Identifier: GNU Lesser General Public Version 2.1
#
# This script formats all the cmake files
#
#
# Usage:
# ./contrib/utilities/cmake_format.sh
#

if test ! -d src -o ! -d include -o ! -d applications; then
	echo "This script must be run from the top-level directory of PRISMS-PF"
	exit 1
fi

if ! [ -x "$(command -v gersemi)" ]; then
	echo "make sure gersemi is in the path"
	exit 1
fi

# Format the root cmake file
gersemi "CMakeLists.txt" -i
echo "Formatted CMakeLists.txt"

# Directories containing the cmake files
TARGET_DIRS=("cmake" "src" "doc" "applications" "tests")

for TARGET_DIR in "${TARGET_DIRS[@]}"; do
	# Check if the directory exists
	if [ ! -d "$TARGET_DIR" ]; then
		echo "Directory $TARGET_DIR does not exist"
		exit 2
	fi

	# Find all cmake files and format excluding the ones 
	# generated in CMakeFiles
	find "$TARGET_DIR" \
		-type d -name CMakeFiles -prune -o \
		-type f \( -name "*.cmake" -o -name "CMakeLists.txt" \) -print | while read -r file; do
		gersemi --in-place "$file"
		echo "Formatted $file"
	done
done

echo "Done"
exit 0
