# Visit 2.13.0 log file
ScriptVersion = "2.13.0"
if ScriptVersion != Version():
    print "This script is for VisIt %s. It may not work with version %s" % (ScriptVersion, Version())
ShowAllWindows()
OpenDatabase("/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/grainGrowth_dream3d/dream3d_test.vtk", 0)
# The UpdateDBPluginInfo RPC is not supported in the VisIt module so it will not be logged.
AddPlot("Pseudocolor", "FeatureIds", 1, 1)
DrawPlots()
