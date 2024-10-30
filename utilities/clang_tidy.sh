#!/bin/bash

#
# Script to run clang-tidy over the entire repo.
#
# Run the script in the main folder (i.e., ./utilities/clang_tidy.sh)
#
# Requires LLVM 16.0.0 or later.

dir=$1

# Check that the user has supplied a directory to begin the search
if [ ! -d "$dir/src" ] || [ ! -d "$dir/include" ] || [ ! -d "$dir/applications" ] || [ ! -f "$dir/CMakeLists.txt" ]; then
    echo "Error, usage is: ./utilities/clang_tidy.sh /path/to/your/library"
    exit 1
fi

# Print source directory
echo "Source directory: ${dir}"

# Check that run-clang-tidy is installed and in path
if ! [ -x "$(command -v run-clang-tidy)" ] || ! [ -x "$(command -v clang++)" ]; then
    echo "Error, run-clang-tidy and clang++ do not seem to be installed/in your path"
    exit 2
fi

# Setup cmake arguments with export compile commands and debug mode
arguments=("-D" "CMAKE_EXPORT_COMPILE_COMMANDS=ON" "-D" "CMAKE_BUILD_TYPE=Debug")

# Compile with clang
CC=clang CXX=clang++ cmake "${arguments[@]}" "$dir" || (echo "cmake failed!"; false) || exit 2
cmake --build . || (echo "make failed!"; false) || exit 3

# Run clang-tidy
run-clang-tidy -p . -quiet -header-filter "$dir/include/(?!IntegrationTools).*" 2>error.txt >output.txt

# Remove duplicates & log them
grep -E '(warning|error): ' output.txt | sort | uniq > clang-tidy.log

# If there are errors, exit unsuccessfully
if [ -s clang-tidy.log ]; then
    cat clang-tidy.log
    exit 4
fi

# Successfully exit
echo "Clang-tidy checks passed!"
exit 0