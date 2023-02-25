#!/usr/bin/env python3

# Script to convert a rectilinear grid vtk file
# from cell data into point data
import sys
from visit import *

#Input file name
ifname="microstructure_RL.vtk"

#Output file name
ofname="microstructure_RL_resampled"

print("Opening input vtk file from DREAM.3D")

# Step 1: Open an Rectilinear vtk file
OpenDatabase(ifname)

# Step 2: Add FeatureIds plot
AddPlot("Pseudocolor", "FeatureIds", 1, 1)
DrawPlots()

#Get the number of grid points in the domain
Query("Grid Information")
gpq = GetQueryOutputValue()
# Extracting number of grid points in each direction
num_x_coords = int(gpq[2])
num_y_coords = int(gpq[3])
num_z_coords = int(gpq[4])

print ("number of grid points z, y, z: " + str(num_x_coords) + ", " + str(num_y_coords) + ", " + str(num_z_coords))

#Getting the min and max coordinates of the system
Query("SpatialExtents", use_actual_data=1)
seq = GetQueryOutputValue()
min_x_val = seq[0]
min_y_val = seq[2]
min_z_val = seq[4]

print("Translation of system coordinades to origin")

# Translating the system such that all coordinates are positive
AddOperator("Transform", 1)
TransformAtts = TransformAttributes()
if min_x_val < 0.0:
    TransformAtts.doTranslate = 1
    TransformAtts.translateX = -min_x_val
if min_y_val < 0.0:
    TransformAtts.doTranslate = 1
    TransformAtts.translateY = -min_y_val
if min_z_val < 0.0:
    TransformAtts.doTranslate = 1
    TransformAtts.translateZ = -min_z_val
SetOperatorOptions(TransformAtts, 0, 1)

print("Resampling into rectilinear grid with point data")
#Resampling the system so that values at points are reported
AddOperator("Resample", 1)
ResampleAtts = ResampleAttributes()
ResampleAtts.samplesX = num_x_coords
ResampleAtts.samplesY = num_y_coords
# If the system is only one cell thick in the 2-direction, convert it to a 2D system
if num_z_coords <= 2:
    ResampleAtts.is3D = 0
else:
    ResampleAtts.is3D = 1
    ResampleAtts.samplesZ = num_z_coords
SetOperatorOptions(ResampleAtts, 1, 1)
# Apply transformations
DrawPlots()

#Export the data in VTK Legacy ASCII format
print("Generating resampled output file (VTK Legacy ASCII; RECTILINEAR_GRID format)")

ExportDBAtts = ExportDBAttributes()
ExportDBAtts.filename = ofname
ExportDBAtts.db_type = "VTK"
ExportDBAtts.opts.types = (5, 0)
ExportDatabase(ExportDBAtts)

#Closing the session
sys.exit()
