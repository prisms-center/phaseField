import sys

variable_1 = 'c'
variable_2 = 'c'

output_num_1 = "1000"
output_num_2 = "1000"

directory_1 = "reference_solution/solution-"+output_num_1+".pvtu"
directory_2 = "run_001/solution-"+output_num_2+".pvtu"


OpenDatabase(directory_1)
DefineScalarExpression("diff","("+variable_1+"-pos_cmfe(<"+directory_2+":"+variable_2+">, mesh,0))")


#AddPlot("Contour","diff")
AddPlot("Pseudocolor","diff")

PseudocolorAtts = PseudocolorAttributes()
PseudocolorAtts.limitsMode = PseudocolorAtts.CurrentPlot
SetPlotOptions(PseudocolorAtts)


DrawPlots()

Query("MinMax", use_actual_data=1)
min_max_val = GetQueryOutputValue()
error = max(abs(min_max_val[0]),abs(min_max_val[1]))

text_file = open("error.txt","w")
text_file.write(str(error)) 
text_file.close() 

sys.exit()