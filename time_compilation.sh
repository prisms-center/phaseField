#!/bin/bash

# Clean up the build directory
make distclean
rm CMakeCache.txt
rm .ninja_log

# Build
cmake . -G Ninja

# Start timing the compilation
echo "Starting compilation..."
start_time=$(date +%s)

# Run the compilation command (e.g., make)
ninja -j4

# End timing the compilation
end_time=$(date +%s)

# Calculate the elapsed time
elapsed_time=$((end_time - start_time))

# Print the elapsed time
echo "Compilation completed in $elapsed_time seconds."

# Run ninjatracing
./../ninjatracing/ninjatracing -a .ninja_log > trace.json