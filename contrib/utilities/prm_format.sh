#!/bin/bash
# SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
# SPDX-License-Identifier: GNU Lesser General Public Version 2.1
#
# This script formats all the .prm files
#
#
# Usage:
# ./contrib/utilities/prm_format.sh
#

if test ! -d src -o ! -d include -o ! -d applications; then
	echo "This script must be run from the top-level directory of PRISMS-PF"
	exit 1
fi

# Directories containing the .prm files
TARGET_DIRS=("tests" "applications")

for TARGET_DIR in "${TARGET_DIRS[@]}"; do
	# Check if the directory exists
	if [ ! -d "$TARGET_DIR" ]; then
		echo "Directory $TARGET_DIR does not exist"
		exit 2
	fi

	# Find all .prm files and format
	find "$TARGET_DIR" -type f \( -name "*.prm" \) | while read -r file; do
		# Run prm_tools.py on the file
		python3 contrib/utilities/prm_tools.py "$file" "$file"
		echo "Formatted $file"
	done
done

echo "Done"
exit 0
