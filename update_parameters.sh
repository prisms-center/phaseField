# !/usr/bin/bash

# Find all parameters.prm files in the applications folder and apply sed
find applications -type f -name "parameters.prm" | while read -r file; do
    sed -i 's/set Number of dimensions = /set dim = /' "$file"
    sed -i 's/set Refine factor = /set global refinement = /' "$file"
    sed -i 's/set Element degree = /set degree = /' "$file"
    echo "Updated $file"
done
