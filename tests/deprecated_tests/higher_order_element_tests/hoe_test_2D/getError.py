import sys

variable_1 = 'c'
variable_2 = 'c'

output_num_1 = "100000"
output_num_2 = "100000"

f = open('run_num_1.txt','r')
run_num_1 = f.read()
f.close()

f = open('run_num_2.txt','r')
run_num_2 = f.read()
f.close()

use_mesh_1 = True
apply_3slice = False
apply_box_crop = False

directory_1 = run_num_1+"/solution-"+output_num_1+".pvtu"
directory_2 = run_num_2+"/solution-"+output_num_2+".pvtu"

if use_mesh_1 == True: 
	OpenDatabase(directory_1)
	DefineScalarExpression("diff","abs(("+variable_1+"-pos_cmfe(<"+directory_2+":"+variable_2+">, mesh,0))^2)")
else:
	OpenDatabase(directory_2)
	DefineScalarExpression("diff","abs(("+variable_2+"-pos_cmfe(<"+directory_1+":"+variable_1+">, mesh,0))^2)")

#AddPlot("Contour","diff")
AddPlot("Pseudocolor","diff")

PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.limitsMode = PseudocolorAtts.CurrentPlot
SetPlotOptions(PseudocolorAtts)

if apply_box_crop == True:
	AddOperator("Box", 1)
	BoxAtts = BoxAttributes()
	BoxAtts.amount = BoxAtts.Some  # Some, All
	BoxAtts.minx = 2.0
	BoxAtts.maxx = 30.0
	BoxAtts.miny = 2.0
	BoxAtts.maxy = 30.0
	BoxAtts.minz = 2.0
	BoxAtts.maxz = 30.0
	BoxAtts.inverse = 0
	SetOperatorOptions(BoxAtts, 1)

if apply_3slice == True:
	AddOperator("ThreeSlice")
	slice_atts = ThreeSliceAttributes()
	slice_atts.x = 16
	slice_atts.y = 16
	slice_atts.z = 16
	SetOperatorOptions(slice_atts)

DrawPlots()

#Query("MinMax", use_actual_data=1)
#min_max_val = GetQueryOutputValue()
#error = max(abs(min_max_val[0]),abs(min_max_val[1]))

# Calculate the L2 error:
Query("Weighted Variable Sum", use_actual_data=1)
error = GetQueryOutputValue()


text_file = open("error.txt","w")
text_file.write(str(error)) 
text_file.close() 

sys.exit()