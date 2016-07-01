import subprocess
import shutil
import glob
import math
import os.path

# ----------------------------------------------------------------------------------------
# Function that generates a header file for the simulation, compiles the PRISMS-PF code
# and runs the executable.
# ----------------------------------------------------------------------------------------
def run_simulation(refineFactor,degree,subdivisionsX,subdivisionsY,subdivisionsZ,run_name):

	text_file = open("space_parameters.h","w")
	text_file.write("// Parameters list for the parameters involved in mesh generation \n") 
	text_file.write("#define subdivisionsX " + str(subdivisionsX) + "\n") 
	text_file.write("#define subdivisionsY " + str(subdivisionsY) + "\n") 
	text_file.write("#define subdivisionsZ " + str(subdivisionsZ) + "\n") 
	text_file.write("#define refineFactor " + str(refineFactor) + "\n") 
	text_file.write("#define finiteElementDegree " + str(degree) + "\n")
	text_file.close()

	subprocess.call(["make", "release"])
	#subprocess.call(["./main"])
	subprocess.call(["mpirun", "-n", "6", "main"])

	subprocess.call(["mkdir",run_name])
	for output_files in glob.glob('*vtu'):
		shutil.move(output_files,run_name)
		
# ----------------------------------------------------------------------------------------
# Function that calculates the maximum error between a run and a reference result.
# ----------------------------------------------------------------------------------------
def getError(run_num):

	text_file = open("run_num_1.txt","w")
	
	if run_num < 10:
		text_file.write('run_00'+str(run_num)) 
	elif run_num < 100:
		text_file.write('run_0'+str(run_num)) 
	elif run_num < 1000:
		text_file.write('run_'+str(run_num))
	else:
		print "error: run number must be lower than 1000"
		
	text_file.close() 	
	subprocess.call(["visit","-nowin", "-cli","-s","getError.py"])
	f = open('error.txt','r')
	error = float(f.read())
	f.close()
	
	return error
		
# ----------------------------------------------------------------------------------------

# Choose whether or not to re-run the reference simulation (warning, takes several hours)
re_run_ref_simulation = False

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
if os.path.exists("run_006") == True:
	shutil.rmtree("run_006")
if re_run_ref_simulation == True:
	if os.path.exists("reference_solution") == True:
		shutil.rmtree("reference_solution")

# Run simulations with varying mesh size and element order:
run_simulation(3,1,1,1,1,'run_001')
run_simulation(4,1,1,1,1,'run_002')
run_simulation(5,1,1,1,1,'run_003')
run_simulation(2,2,1,1,1,'run_004')
run_simulation(3,2,1,1,1,'run_005')
run_simulation(4,2,1,1,1,'run_006')
if re_run_ref_simulation == True:
	run_simulation(7,1,1,1,1,'reference_solution')


# Find the error for runs 1-6 (as compared to the reference solution)
text_file = open("run_num_2.txt","w")
text_file.write('reference_solution') 
text_file.close() 	

error = []
error.append(getError(1))
error.append(getError(2))
error.append(getError(3))
error.append(getError(4))
error.append(getError(5))
error.append(getError(6))

# Calculate and output the error

print "Error results: \n"
print "First order element error:", str(error[0]), "(refine level 3)", str(error[1]), "(refine level 4)", \
	str(error[2]), "(refine level 5)", "\n"

print "Second order element error:", str(error[3]), "(refine level 2)", str(error[4]), "(refine level 3)", \
	str(error[5]), "(refine level 4)", "\n"
