#!/bin/bash

#
# This script updates the copyright information in all .cc and .h files based on 
# the year that file was lasted updated (according to git).
#
#
# Usage:
# ./contrib/utilities/update_copyright.sh
#

if test ! -d src -o ! -d include -o ! -d applications ; then
  echo "This script must be run from the top-level directory of PRISMS-PF"
  exit 0
fi

# Directories containing the .cc and .h files
TARGET_DIRS=("src" "include" "tests" "applications")

for TARGET_DIR in "${TARGET_DIRS[@]}"; do
  # Check if the directory exists
  if [ ! -d "$TARGET_DIR" ]; then
    echo "Directory $TARGET_DIR does not exist."
    exit 1
  fi

  # Find all .cc and .h files and insert the text
  find "$TARGET_DIR" -type f \( -name "*.cc" -o -name "*.h" \) | while read -r file; do
    last_year=`git log -n 3 --date=short --format="format:%cd %s" ${file} | \
      egrep -i -v "update.*copyright|copyright.*update|Update license headers" | \
      head -n 1 | \
      perl -p -e 's/^(\d\d\d\d)-.*/\1/g;'`
    INSERT_TEXT="// SPDX-FileCopyrightText: © ${last_year} PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1\n\n"
    # Insert the text at the beginning of the file
    echo -e "$INSERT_TEXT$(cat "$file")" > "$file"
    echo "Inserted text into $file"
  done
done

echo "Done"