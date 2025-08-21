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

FILE_PATH="$1"

# Grab the latest version from the VERSION file
latest_version=$(cat VERSION | tr -d ' \n')
echo "Latest version found: $latest_version"

# Create HTML
cat << EOF > $FILE_PATH/index.html
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta http-equiv="refresh" content="0; url=$latest_version/index.html">
    <title>Redirecting...</title>
</head>
<body>
    <p>If you are not redirected automatically, <a href="$latest_version/index.html">click here</a>.</p>
</body>
</html>
EOF

echo "Done"
exit 0
