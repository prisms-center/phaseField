#!/bin/bash
# SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
# SPDX-License-Identifier: GNU Lesser General Public Version 2.1
#
# This script updates the copyright information in all .cc and .h files based on
# the year that file was lasted updated (according to git).
#
#
# Usage:
# ./contrib/utilities/update_copyright.sh
#

if test ! -d src -o ! -d include -o ! -d applications; then
  echo "This script must be run from the top-level directory of PRISMS-PF"
  exit 1
fi

if command -v clang-format &>/dev/null; then
  CLANG_FORMAT=clang-format
elif command -v clang-format-18 &>/dev/null; then
  CLANG_FORMAT=clang-format-18
else
  echo "Neither clang-format nor clang-format-18 are in path."
  exit 1
fi

# Directories containing the .cc and .h files
TARGET_DIRS=("src" "include" "tests" "applications")

for TARGET_DIR in "${TARGET_DIRS[@]}"; do
  # Check if the directory exists
  if [ ! -d "$TARGET_DIR" ]; then
    echo "Directory $TARGET_DIR does not exist"
    exit 2
  fi

  # Find all .cc and .h files and insert the text
  find "$TARGET_DIR" -type f \( -name "*.cc" -o -name "*.h" \) | while read -r file; do
    last_year=$(git log -n 3 --date=short --format="format:%cd %s" ${file} |
      egrep -i -v "update.*copyright|copyright.*update|Update license headers" |
      head -n 1 |
      perl -p -e 's/^(\d\d\d\d)-.*/\1/g;')
    INSERT_TEXT="// SPDX-FileCopyrightText: © ${last_year} PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1
"
    # Check if the file already contains the copyright text
    if grep -q "SPDX-FileCopyrightText" "$file"; then
      # Update the year in the existing text
      sed -i "s|// SPDX-FileCopyrightText: © .* PRISMS Center at the University of Michigan|// SPDX-FileCopyrightText: © ${last_year} PRISMS Center at the University of Michigan|g" "$file"
      echo "Updated text in $file"
    else
      # Insert the text at the beginning of the file
      printf "%s\n%s" "$INSERT_TEXT" "$(cat "$file")" >"$file"
      echo "Inserted text into $file"
    fi

    # Run clang-format on the file
    $CLANG_FORMAT -i "$file"
    echo "Formatted $file with clang-format"
  done
done

echo "Done"
exit 0
