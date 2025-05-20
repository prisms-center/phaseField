#!/bin/bash
# SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
# SPDX-License-Identifier: GNU Lesser General Public Version 2.1
#
# This script removes all comments from all .prm files
#
#
# Usage:
# ./contrib/utilities/prm_remove_comments.sh
#

if test ! -d src -o ! -d include -o ! -d applications; then
	echo "This script must be run from the top-level directory of PRISMS-PF"
	exit 1
fi

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
MAIN_DIR="$SCRIPT_DIR/../../"

find "$MAIN_DIR" -type f -name "*.prm" -exec sed -i -E 's/^\s*#.*//; s/\s+#.*//' {} +

echo "Done"
exit 0
