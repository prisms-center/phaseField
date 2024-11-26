#!/usr/bin/env python3

import sys
from visit import *

#DeleteAllPlots()

# Step 1: Open a database (the whole .vtu time series)
dbname="solution-*.vtu database"
OpenDatabase(dbname)

#Step 2: Add contour plot to obtain grid information and dimensionality
#Add Contour plot
AddPlot("Contour", "n")
ContourAtts = ContourAttributes()
ContourAtts.contourValue = (0.5)
ContourAtts.contourMethod = ContourAtts.Value
SetPlotOptions(ContourAtts)
# Step 3: Draw the plots
DrawPlots()
        
# Step 4: Get the spacial extents
Query("SpatialExtents")
gpq = GetQueryOutputValue()
# Extracting number of coordinate limits
#(i.e. min x, max x, min y, max y, min z, max z
dimslength=len(gpq);
        
# Step 5 Determine if system is 2D or 3D (2D - dimslength=4; 3D - dimslength=6)
if dimslength >= 5:
    is2D = False
else:
    is2D = True
        
if is2D:
    DeleteActivePlots()
    AddPlot("Pseudocolor", "n")
    DrawPlots()

# Step 7: Animate through time and save results
for states in range(TimeSliderGetNStates()):
    #Set slider to state
    SetTimeSliderState(states)
    # Get the time corresponding to the state
    Query("Time")
    # Assign this time to the variable "t"
    t = GetQueryOutputValue()
    # Print the state number, time and phase fraction to
    # screen and to files
    print("Saving frame % d, time %.1f" %(states, t))
    SaveWindowAtts = SaveWindowAttributes()
    SaveWindowAtts.fileName = "frame_"
    SetSaveWindowAttributes(SaveWindowAtts)
    SaveWindow()
DeleteAllPlots()
CloseDatabase(dbname)

sys.exit()
