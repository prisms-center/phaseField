#!/usr/bin/env python3

import sys

from numpy import *
from visit import *

# DeleteAllPlots()

# Step 1: Open a database (the whole .vtu time series)
dbname = "solution-*.vtu database"
OpenDatabase(dbname)

# Defining number of samples in each direction (for resampling operator)
# The number of samples in the direction in which the system is the largest
# gets assigned nmax samples. The remaining directions get assigned
# a propostional number according to their length
# Thus number can be increased if domains are underresolved or
# decreased if domains are sufficiently resolved and calculations
# take too long
nmax = 256

# Step 2: Add contour plot to obtain grid information and dimensionality
# Add Contour plot
AddPlot("Contour", "n")
ContourAtts = ContourAttributes()
ContourAtts.contourValue = 0.5
ContourAtts.contourMethod = ContourAtts.Value
SetPlotOptions(ContourAtts)
# Step 3: Draw the plots
DrawPlots()

# Step 4: Get the spacial extents
Query("SpatialExtents")
gpq = GetQueryOutputValue()
# Extracting number of coordinate limits
# (i.e. min x, max x, min y, max y, min z, max z
dimslength = len(gpq)

# Step 5 Determine if system is 2D or 3D (2D - dimslength=4; 3D - dimslength=6)
if dimslength >= 5:
    is3D = True
else:
    is3D = False

DeleteActivePlots()

# Preparing system for resampling
AddPlot("Pseudocolor", "n")
DrawPlots()

# Step 3: Get the number of grid points in the domain
Query("Grid Information")
gpq = GetQueryOutputValue()
# Extracting number of grid points in the z- direction
num_z_coords = int(gpq[4])

# Step 3: Get the system length in each direction
Query("SpatialExtents")
sq = GetQueryOutputValue()
lx = sq[1] - sq[0]
ly = sq[3] - sq[2]
if is3D == True:
    lz = sq[5] - sq[4]
else:
    lz = 0.0

lv = [lx, ly, lz]

# Step 4: Calculate the numbers of points in each direction for the
# resample operators
nsx = int(nmax * lv[0] / max(lv))
nsy = int(nmax * lv[1] / max(lv))
nsz = int(nmax * lv[2] / max(lv))

# Step 5: Resample the data to create a strucutured mesh
# Adding the Resample Operator
AddOperator("Resample")
# The object r contains the Resample Operator Attributes
r = ResampleAttributes()
# Setting the number of samples in x and y
# (this may need to be adjusted for different system sizes
# and refine factors)
r.samplesX = nsx
r.samplesY = nsy
# For 3D r.samplesZ can be added
if is3D == True:
    r.is3D = 1
    r.samplesZ = nsz
else:
    r.is3D = 0

# Setting the operator options
SetOperatorOptions(r)

# Step 6: Split order parameter domains into different regions
# based on value
# Adding operator Isovolume
AddOperator("Isovolume")
# The object iv contains the Isovolume Operator Attributes
iv = IsovolumeAttributes()
# Setting the lower bound to 0.5
# This will select all the regions for which n>0.5
iv.lbound = 0.5
# Setting the operator options
SetOperatorOptions(iv)
DrawPlots()

# Step 7: Initializing the number of domains for each time state
# And opening output file
nodoms = [0.0] * TimeSliderGetNStates()
# Setting the output file name
ofnm = "domain_stats_vs_t.txt"
# opening file to output
outF = open(ofnm, "w")
print >> outF, "# frame no. / time / no. domains / avg. domain vol. (area) / std. dev."

print("# frame no. / time / no. domains / avg. domain vol. (area) / std. dev.")

# Step 8: Animate through time and save results
# Count number of domains and save output
for states in range(TimeSliderGetNStates()):
    # Set slider to state
    SetTimeSliderState(states)
    # Get the time corresponding to the state
    Query("Time")
    # Assign this time to the variable "t"
    t = GetQueryOutputValue()
    # Get the number of domains as the number of connected
    # components
    Query("Number of Connected Components")
    # Assign the output of the query to the variable nodoms
    # for the current state
    nodoms = GetQueryOutputValue()
    if nodoms > 0:
        Query("Connected Components Summary")
        res = GetQueryOutputObject()
        # Show names in the results dictionary
        # print(res.keys())
        # calculate the volume (or area) of each connected component
        if r.is3D:
            # calculate the area of each connected component
            size_vec = res["comp_volume"]
        else:
            size_vec = res["comp_area"]
        # Convert tuple to numpy array
        size_vec_np = asarray(size_vec)
        # Average domains size
        avg_size = average(size_vec_np)
        # Domains size std. deviation
        std_dev = std(size_vec_np)
    else:
        avg_size = 0.0
        std_dev = 0.0
    # Print the state number, time and no. of domains to
    # screen and to files
    print("% d, %.7f, % d, %.7f, %.7f" % (states, t, nodoms, avg_size, std_dev))
    print >> outF, "% d, %.7f, % d, %.7f, %.7f" % (states, t, nodoms, avg_size, std_dev)
# Close output file
outF.close()
# Delete all plots
DeleteAllPlots()
# Close the database
CloseDatabase(dbname)

sys.exit()
