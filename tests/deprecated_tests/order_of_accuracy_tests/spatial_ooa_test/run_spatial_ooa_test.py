import subprocess
import shutil
import glob
import math
import os.path

# ----------------------------------------------------------------------------------------
# Function that generates a header file for the simulation, compiles the PRISMS-PF code
# and runs the executable.
# ----------------------------------------------------------------------------------------
def run_simulation(refineFactor,run_name):

	subdivisionsX = 1
	subdivisionsY = 1
	subdivisionsZ = 1
	finiteElementDegree = 1

	text_file = open("space_parameters.h","w")
	text_file.write("// Parameters list for the parameters involved in mesh generation \n") 
	text_file.write("#define subdivisionsX " + str(subdivisionsX) + "\n") 
	text_file.write("#define subdivisionsY " + str(subdivisionsY) + "\n") 
	text_file.write("#define subdivisionsZ " + str(subdivisionsZ) + "\n") 
	text_file.write("#define refineFactor " + str(refineFactor) + "\n") 
	text_file.write("#define finiteElementDegree " + str(finiteElementDegree) + "\n")
	text_file.close()

	subprocess.call(["make", "release"])
	#subprocess.call(["./main"])
	subprocess.call(["mpirun", "-n", "4", "main"])

	subprocess.call(["mkdir",run_name])
	for output_files in glob.glob('*vtu'):
		shutil.move(output_files,run_name)
		
# ----------------------------------------------------------------------------------------

# If files exist from previous tests, delete them
if os.path.exists("run_001") == True:
	shutil.rmtree("run_001")
if os.path.exists("run_002") == True:
	shutil.rmtree("run_002")
if os.path.exists("run_003") == True:
	shutil.rmtree("run_003")
if os.path.exists("run_004") == True:
	shutil.rmtree("run_004")

# Run simulations with a decreasing time step (and an corresponding increase in number of iterations)
run_simulation(4,'run_001')
run_simulation(5,'run_002')
run_simulation(6,'run_003')
run_simulation(7,'run_004')

# Find the error for runs 1-3 (as compared to run 4)	
subprocess.call(["visit", "-cli","-s","getError1.py"])
subprocess.call(["visit", "-cli","-s","getError2.py"])
subprocess.call(["visit", "-cli","-s","getError3.py"])

f = open('error1.txt','r')
error1 = float(f.read())
f.close()

f = open('error2.txt','r')
error2 = float(f.read())
f.close()

f = open('error3.txt','r')
error3 = float(f.read())
f.close()

# Calculate and output order of accuracy
ooa_1_to_2 = math.log10(error1/error2)/math.log10(2)
ooa_2_to_3 = math.log10(error2/error3)/math.log10(2)

print "Order of accuracy results: \n"
print "Mesh refinement of 4 to mesh refinement of 5:", str(ooa_1_to_2), "\n"
print "Mesh refinement of 5 to mesh refinement of 6:", str(ooa_2_to_3), "\n"

text_file = open("spat_ooa_results.txt","w")
text_file.write("Results of standard spatial order of accuracy test: \n") 
text_file.write("From coarse mesh to medium mesh: " + str(ooa_1_to_2) + "\n") 
text_file.write("From from medium mesh to fine mesh: " + str(ooa_2_to_3) + "\n") 
text_file.close()

