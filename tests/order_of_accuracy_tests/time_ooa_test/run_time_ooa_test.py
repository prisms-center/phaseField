import subprocess
import shutil
import glob
import math
import os.path

# ----------------------------------------------------------------------------------------
# Function that generates a header file for the simulation, compiles the PRISMS-PF code
# and runs the executable.
# ----------------------------------------------------------------------------------------
def run_simulation(factor,run_name):

	timeIncrements = int(100*factor)
	timeFinal = 100000000
	skipImplicitSolves = 10000000
	timeStep = 1.5e-4/factor

	text_file = open("time_parameters.h","w")
	text_file.write("// Parameters list for the parameters involved in time stepping \n") 
	text_file.write("#define timeStep " + str(timeStep) + "\n") 
	text_file.write("#define timeIncrements " + str(timeIncrements) + "\n") 
	text_file.write("#define timeFinal " + str(timeFinal) + "\n") 
	text_file.write("#define skipImplicitSolves " + str(skipImplicitSolves) + "\n") 
	text_file.close()

	subprocess.call(["make", "release"])
	subprocess.call(["./main"])

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
run_simulation(1.0,'run_001')
run_simulation(2.0,'run_002')
run_simulation(4.0,'run_003')
run_simulation(16.0,'run_004')

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
print "Time step 1 to time step 2:", str(ooa_1_to_2), "\n"
print "Time step 2 to time step 3:", str(ooa_2_to_3), "\n"

text_file = open("time_ooa_results.txt","w")
text_file.write("Results of standard temporal order of accuracy test: \n") 
text_file.write("From time step 1 to time step 2: " + str(ooa_1_to_2) + "\n") 
text_file.write("From time step 2 to time step 3: " + str(ooa_2_to_3) + "\n") 
text_file.close()

