#!/bin/bash
# SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
# SPDX-License-Identifier: GNU Lesser General Public Version 2.1
#
# This script unwraps the compile_commands.json for the mpicxx compiler wrapper.
# The compile wrapper can cause issues with clang-tidy as described here:
# https://github.com/clangd/clangd/issues/2328
#
# Instead of using the mpicxx wrapper, we use the actual compiler (e.g., g++) and
# all passed arguments as you see with mpicxx -show. For now this script is only
# enable when a specific CMake flag is passed (i.e., UNWRAP_COMPILER=ON).
#
#
# Usage:
# ./contrib/utilities/unwrap_compile_commands.sh
#

if test ! -d src -o ! -d include -o ! -d applications; then
    echo "This script must be run from the top-level directory of PRISMS-PF"
    exit 1
fi

# Check that compile_commands.json exists
if test ! -f compile_commands.json; then
    echo "There is no compile_commands.json file in the current directory."
    exit 1
fi

# Make a copy of compile commands for us to operate on
cp compile_commands.json compile_commands.json.bak || {
    echo "Failed to copy compile_commands.json"
    exit 1
}

# Grab the compiler and unwrap it
COMPILER=$(grep -m1 '"command":' compile_commands.json | sed -E 's/.*"command": ?"([^"]+)".*/\1/' | awk '{print $1}') || {
    echo "Compiler not found or not executable: $COMPILER"
    exit 1
}
UNWRAPPED_COMPILER=$($COMPILER -show) || {
    echo "Failed to run '$COMPILER -show'"
    rm compile_commands.json.bak
    exit 0
}

# Replace the COMPILER instances with the UNWRAPPED_COMPILER
sed -i "s|$COMPILER|$UNWRAPPED_COMPILER|g" compile_commands.json || {
    echo "sed failed to update compile_commands.json"
    exit 1
}

# Delete the backup file
rm compile_commands.json.bak || {
    echo "Failed to delete backup"
    exit 1
}

exit 0
