#!/usr/bin/env python3

import sys
from visit import *
from numpy import *

#DeleteAllPlots()

#Defining number of samples in each direction (for resampling operator)
#The number of samples in the direction in which the system is the largest
#gets assigned nmax samples. The remaining directions get assigned
#a propostional number according to their length
#Thus number can be increased if domains are underresolved or
#decreased if domains are sufficiently resolved and calculations
#take too long
nmax=1024;

# Step 1: Open a database (the whole .vtu time series)
dbname="solution-*.vtu database"
OpenDatabase(dbname)

# Step 2: Add pseudocolor plot (using variable "n")
# This variable must be in the range [0,1]
# with 0 representing one phase, 1 representing another phase
# and n=0.5 representing the midpoint accross the interface
AddPlot("Pseudocolor", "n")
DrawPlots()

# Step 3: Get the number of grid points in the domain
Query("Grid Information")
gpq = GetQueryOutputValue()
# Extracting number of grid points in the z- direction
num_z_coords = int(gpq[4])

#Step 3: Get the system length in each direction
Query("SpatialExtents")
sq=GetQueryOutputValue()
lx=sq[1]-sq[0]
ly=sq[3]-sq[2]
if num_z_coords >= 1:
    lz=sq[5]-sq[4]
else:
    lz = 0.0
lv=[lx,ly,lz]

#Step 4: Calculate the numbers of points in each direction for the
#resample operators
nsx=int(nmax*lv[0]/max(lv))
nsy=int(nmax*lv[1]/max(lv))
nsz=int(nmax*lv[2]/max(lv))

# Step 5: Resample the data to create a strucutured mesh
# Adding the Resample Operator
AddOperator("Resample")
#The object r contains the Resample Operator Attributes
r=ResampleAttributes()
# Setting the number of samples in x and y
# (this may need to be adjusted for different system sizes
# and refine factors)
r.samplesX = nsx
r.samplesY = nsy
# For 3D r.samplesZ can be added
if num_z_coords >= 1:
    r.is3D = 1
    r.samplesZ = nsz
else:
    r.is3D=0

# Setting the operator options
SetOperatorOptions(r)

# Step 6: Split order parameter domains into different regions
# based on value
# Adding operator Isovolume
AddOperator("Isovolume")
#The object iv contains the Isovolume Operator Attributes
iv=IsovolumeAttributes()
# Setting the lower bound to 0.5
# This will select all the regions for which n>0.5
iv.lbound=0.5
# Setting the operator options
SetOperatorOptions(iv)
DrawPlots()

#Step 7: Initializing the number of domains for each time state
# And opening output file
nodoms=[0.0]*TimeSliderGetNStates()
#Setting the output file name
ofnm="nodoms_vs_t.txt"
#opening file to output
outF = open(ofnm, "w")

# Step 8: Animate through time and save results
# Count number of domains and save output
for states in range(TimeSliderGetNStates()):
    #Set slider to state
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
    nodoms[states]=GetQueryOutputValue()
    # Print the state number, time and no. of domains to
    # screen and to files
    print("% d, %.1f, % d" %(states, t, nodoms[states]))
    print >> outF, "% d %.1f % d" %(states, t, nodoms[states])
# Close output file
outF.close()
# Delete all plots
DeleteAllPlots()
# Close the database
CloseDatabase(dbname)

sys.exit()
