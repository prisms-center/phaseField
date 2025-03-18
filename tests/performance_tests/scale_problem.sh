#!/bin/bash

#
# This script copies a given equation n times in the performance tests.
#
#
# Usage:
# ./scale_problem.sh /APP_DIR N_COPIES
#    with:
#      APP_DIR pointing toward the application directory to update
#      N_COPIES the number of copies of the equation that should be made
#

# Grab the inputs for APP_DIR and N_COPIES
APP_DIR=$1
APP_DIR=$(cd "$APP_DIR";pwd)
N_COPIES=$2

# Check that the paths are correct
if [ ! -f "$APP_DIR/custom_pde.h" ] || [ ! -f "$APP_DIR/CMakeLists.txt" ] ; then
    echo "Usage:"
    echo "  scale_problem.sh /path/to/application <number of copies>"
    exit 1
fi
echo "APP-DIR=$APP_DIR"
echo "N COPIES=$N_COPIES"

# Check that the number of copies is greater than 1
if [ "$N_COPIES" -lt 1 ] ; then
    echo "Error: N_COPIES must be greater than or equal to 1"
    exit 2
fi

# Update the n_copies value in custom_pde.h
sed -i "s|const unsigned int n_copies\s*=\s*[0-9]*\s*;|const unsigned int n_copies = $N_COPIES;|" "$APP_DIR/custom_pde.h"

# Update the parameters file
grep -E "set boundary condition for [a-z]+0\s*=\s*[A-Z]+" "$APP_DIR/parameters.prm" > "$APP_DIR/parameters_new.prm"
sed -i -E "/set boundary condition for [a-z]+[0-9]+\s*=\s*[A-Z]+/d" "$APP_DIR/parameters.prm"
for ((i=0; i<N_COPIES; i++)) ; do
    sed "s/\([a-z]\)0/\1$i/" "$APP_DIR/parameters_new.prm" >> "$APP_DIR/parameters.prm"
done
rm "$APP_DIR/parameters_new.prm"

# Compile and run 
cd "$APP_DIR"
make -j$(nproc)
for ((i=0; i<3; i++)) ; do
    mpirun -n 1 main -P runtime-report,mem.highwatermark > "trial_${i}_${N_COPIES}.txt" 2>&1
done