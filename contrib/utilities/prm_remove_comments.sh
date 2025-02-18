#!/bin/bash

#
# This script removes all comments from all .prm files
#
#
# Usage:
# ./contrib/utilities/prm_remove_comments.sh
#

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
MAIN_DIR="$SCRIPT_DIR/../../"

find "$MAIN_DIR" -type f -name "*.prm" -exec sed -i -E 's/^\s*#.*//; s/\s+#.*//' {} +