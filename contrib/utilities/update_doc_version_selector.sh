#!/bin/bash
# SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
# SPDX-License-Identifier: GNU Lesser General Public Version 2.1
#
# This script updates the version selector in the documentation.
#
#
# Usage:
# ./contrib/utilities/update_doc_version_selector.sh FILE_PATH
#

REMOTE="origin"
BRANCH="gh-pages"
FILE_PATH="$1"

# Find the files in the doxygen folder of the remote branch
git fetch "$REMOTE" "$BRANCH"
dirs=$(git ls-tree --name-only -d "$REMOTE/$BRANCH":doxygen | sort -rV)
echo "Found directories under doxygen/:"
echo "$dirs"
echo ""

# Get the current git branch
current_branch=$(git rev-parse --abbrev-ref HEAD)
echo "Current branch: $current_branch"

# Add the branch to the top of the list if it's not already there
if ! echo "$dirs" | grep -q "^$current_branch$"; then
  dirs="$current_branch"$'\n'"$dirs"
fi
echo "Final list of versions:"
echo "$dirs"
echo ""

# Create HTML
echo '<select id="versionSelector">' >$FILE_PATH/version_selector.html
for dir in $dirs; do
  if [[ "$(basename "$dir")" != .* ]]; then
    version=$(basename "$dir")
    echo "    <option value=\"$version\">$version</option>" >>$FILE_PATH/version_selector.html
  fi
done
echo '</select>' >>$FILE_PATH/version_selector.html

echo "Done"
exit 0
