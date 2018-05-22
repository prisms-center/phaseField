directory = "/Users/stephendewitt/Documents/workspace/PRISMS_workspace/phaseField_dev/applications/grainGrowth_dream3d/"

filename = directory + "dream3d_test.vtk"

OpenDatabase(filename)

AddPlot("Pseudocolor","FeatureIds")
DrawPlots()

opts = GetExportOptions("VTK")
print opts
