#!/bin/bash

# Directory containing the .cc and .h files
TARGET_DIR="$1"

# Check if the directory is provided
if [ -z "$TARGET_DIR" ]; then
  echo "Usage: $0 <directory>"
  exit 1
fi

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
  INSERT_TEXT="// SPDX-FileCopyrightText: © ${last_year} PRISMS Center at the University of Michigan <prisms-pf@umich.edu>
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1\n\n"
  # Insert the text at the beginning of the file
  echo -e "$INSERT_TEXT$(cat "$file")" > "$file"
  echo "Inserted text into $file"
done

echo "Text insertion complete."