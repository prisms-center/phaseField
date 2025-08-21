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

# Grab the latest version from the VERSION file
latest_version=$(cat VERSION | tr -d ' \n')
echo "Latest version found: $latest_version"

# Add the latest version to the top of the list if it's not already there
if ! echo "$dirs" | grep -q "^$latest_version$"; then
	dirs="$latest_version"$'\n'"$dirs"
fi
echo "Final list of versions:"
echo "$dirs"
echo ""

# Create HTML
echo '<select id="versionSelector">' > $FILE_PATH/version_selector.html
for dir in $dirs; do
	if [[ "$(basename "$dir")" != .* ]]; then
		version=$(basename "$dir")
		echo "    <option value=\"$version\">$version</option>" >> $FILE_PATH/version_selector.html
	fi
done
echo '</select>' >> $FILE_PATH/version_selector.html

echo "Done"
exit 0
