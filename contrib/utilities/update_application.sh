#!/bin/bash

#
# This script updates a PRISMS-PF application from an
# older version to the current version.
#
#
# Usage:
# ./contrib/utilities/update_application.sh SRC_DIR /applications/APP_DIR
#    with:
#      SRC_DIR pointing toward the PRISMS-PF source directory
#      APP_DIR pointing toward the application directory to update
#

# Grab the inputs for SRC_DIR and APP_DIR
SRC_DIR=$1
SRC_DIR=$(cd "$SRC_DIR";pwd)
APP_DIR=$2
APP_DIR=$(cd "$APP_DIR";pwd)

# Check that the paths are correct
if test ! -f "$APP_DIR/customPDE.h" -o ! -f "$APP_DIR/../main.cc" -o ! -d "$SRC_DIR/src" -o ! -d "$SRC_DIR/include" -o ! -d "$SRC_DIR/applications" -o ! -f "$SRC_DIR/CMakeLists.txt" ; then
    echo "Usage:"
    echo "  update_application.sh /path/to/PRISMS-PF /path/to/application"
    exit 1
fi
echo "SRC-DIR=$SRC_DIR"
echo "APP-DIR=$APP_DIR"

# Check that the template application (allenCahn) is there
TEMPLATE_DIR=$(realpath "$APP_DIR/../allenCahn")

if test ! -d "$TEMPLATE_DIR" ; then
    echo ""
    echo "The allenCahn application does not exist."
    echo "This is needed to template the update."
    exit 2
fi
echo "TEMPLATE-DIR=$TEMPLATE_DIR"

# Copy CMakeLists.txt
cp "$TEMPLATE_DIR/CMakeLists.txt" "$APP_DIR"

# Remove main.cc
rm "$APP_DIR/main.cc"

# Edit customPDE.h
sed -i 's|dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc|const Point<dim, VectorizedArray<double>> q_point_loc, const VectorizedArray<double> element_volume|g' "$APP_DIR/customPDE.h"

sed -i 's|const const| const|g' "$APP_DIR/customPDE.h"

# Edit equations.cc
sed -i 's|dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc|const Point<dim, VectorizedArray<double>> q_point_loc, const VectorizedArray<double> element_volume|g' "$APP_DIR/equations.cc"

# Edit postprocess.cc, if it exists
if test -f "$APP_DIR/postprocess.cc" ; then
    sed -i 's|dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc|const Point<dim, VectorizedArray<double>> q_point_loc, const VectorizedArray<double> element_volume|g' "$APP_DIR/postprocess.cc"

    sed -i 's|const const| const|g' "$APP_DIR/postprocess.cc"
fi


# Try compiling the application
if cd "$APP_DIR" ; then
    echo "Navigating to $APP_DIR"
fi
if cmake . ; then
    echo "CMake configuration successful"
else
    echo "CMake configuration failed"
    exit 3
fi
if make ; then
    echo "Build successful"
else
    echo "Build failed"
    exit 4
fi


echo "application successfully converted"
exit 0