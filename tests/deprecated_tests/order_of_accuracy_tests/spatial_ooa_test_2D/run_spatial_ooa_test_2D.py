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
# Function that calculates the maximum error between a run and a reference result.
# ----------------------------------------------------------------------------------------
def getError(run_num):

	text_file = open("run_num_1.txt","w")
	text_file.write('run_00'+str(run_num)) 
	text_file.close() 	
	subprocess.call(["visit", "-cli","-s","getError.py"])
	f = open('error.txt','r')
	error = float(f.read())
	f.close()
	
	return error
		
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
if os.path.exists("run_005") == True:
	shutil.rmtree("run_005")

# Run simulations with a decreasing time step (and an corresponding increase in number of iterations)
run_simulation(4,'run_001')
run_simulation(5,'run_002')
run_simulation(6,'run_003')
run_simulation(7,'run_004')
run_simulation(8,'run_005')

# Find the error for runs 1-4 (as compared to run 5)
text_file = open("run_num_2.txt","w")
text_file.write('run_005') 
text_file.close() 	

error = []
error.append(getError(1))
error.append(getError(2))
error.append(getError(3))
error.append(getError(4))

# Calculate and output order of accuracy
ooa_1_to_2 = math.log10(error[0]/error[1])/math.log10(2)
ooa_2_to_3 = math.log10(error[1]/error[2])/math.log10(2)
ooa_3_to_4 = math.log10(error[2]/error[3])/math.log10(2)

print "Order of accuracy results: \n"
print "Mesh refinement of 4 to mesh refinement of 5:", str(ooa_1_to_2), "\n"
print "Mesh refinement of 5 to mesh refinement of 6:", str(ooa_2_to_3), "\n"
print "Mesh refinement of 6 to mesh refinement of 7:", str(ooa_3_to_4), "\n"

text_file = open("spat_ooa_results.txt","w")
text_file.write("Results of 2D spatial order of accuracy test: \n") 
text_file.write("From coarse mesh to medium mesh: " + str(ooa_1_to_2) + "\n") 
text_file.write("From from medium mesh to fine mesh: " + str(ooa_2_to_3) + "\n") 
text_file.write("From from fine mesh to very fine mesh: " + str(ooa_3_to_4) + "\n") 
text_file.close()

