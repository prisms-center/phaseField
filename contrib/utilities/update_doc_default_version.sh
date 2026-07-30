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
BASE_BRANCH="master"
DOC_BRANCH="gh-pages"
FILE_PATH="$1"

# Find the files in the doxygen folder of the remote branch
git fetch "$REMOTE" "$DOC_BRANCH"
dirs=$(git ls-tree --name-only -d "$REMOTE/$DOC_BRANCH":doxygen | sort -rV)
echo "Found directories under doxygen/:"
echo "$dirs"
echo ""

# Get the most recent git tag
latest_tag=$(git describe --tags --abbrev=0 "$REMOTE/$BASE_BRANCH")
latest_tag=${latest_tag#v}
echo "Latest tag: $latest_tag"

# Check that the folder exists on the remote branch (gh-pages)
if ! echo "$dirs" | grep -Fxq "${latest_tag}"; then
  echo "Error: '${latest_tag}' does not exist under doxygen/ on $REMOTE/$DOC_BRANCH"
  exit 1
fi

# Create HTML
cat <<EOF >$FILE_PATH/index.html
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta http-equiv="refresh" content="0; url=$latest_tag/index.html">
    <title>Redirecting...</title>
</head>
<body>
    <p>If you are not redirected automatically, <a href="$latest_tag/index.html">click here</a>.</p>
</body>
</html>
EOF

echo "Done"
exit 0
