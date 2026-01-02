#!/usr/bin/env python3

# Calculates the total area of the interface

# For 2D systems, it calculates the total intertface length

import sys

from visit import *

# Step 1: Open a database (all time steps)
db = "solution-*.vtu database"
OpenDatabase(db)

# Step 2: Add Contour plot (using variable "n")
# This variable must be in the range [0,1]
# with 0 representing one phase, 1 representing another phase
# and n=0.5 representing the midpoint accross the interface
AddPlot("Contour", "n", 1, 1)
ContourAtts = ContourAttributes()
ContourAtts.contourValue = 0.5
ContourAtts.contourMethod = ContourAtts.Value
SetPlotOptions(ContourAtts)

# Step 3: Draw the plot
DrawPlots()

# Step 5: Animate through time and save images
intarea = [0.0] * TimeSliderGetNStates()
ofnm = "iarea_vs_t.txt"
outF = open(ofnm, "w")

for states in range(TimeSliderGetNStates()):
    SetTimeSliderState(states)
    Query("Time")
    t = GetQueryOutputValue()
    qresp = Query("Weighted Variable Sum")
    wvs = GetQueryOutputValue()
    # Weighted Variable Sum integrates the value of the order
    # parameter along the area (length) of the contour
    # surface (curve).
    # Since the countour is taken at value n=0.5,
    # We multiply by 2 to obtain the area (or length)
    if wvs == ():
        intarea[states] = 0.0
    else:
        intarea[states] = 2.0 * wvs
    print("% d, %.1f, %.5f" % (states, t, intarea[states]))
    print >> outF, "% d %.1f %.5f" % (states, t, intarea[states])
outF.close()
DeleteAllPlots()
CloseDatabase(db)

sys.exit()
