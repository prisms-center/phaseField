import subprocess
import shutil
import glob
import math
import os.path
import os

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
if os.path.exists("simulation_result") == True:
	shutil.rmtree("simulation_result")
if os.path.exists("analytical_result") == True:
	shutil.rmtree("analytical_result")

# Run the 1D diffusion simulation
#with cd("diffusionAnalyticComparison"):
os.chdir("diffusionAnalyticComparison")
subprocess.call(["rm", "CMakeCache.txt"])
subprocess.call(["cmake", "."])
subprocess.call(["make", "debug"])
subprocess.call(["./main"])

# Get the analytic result
#with cd("../diffusionAnalyticResult"):
os.chdir("../diffusionAnalyticResult")
subprocess.call(["rm", "CMakeCache.txt"])
subprocess.call(["cmake", "."])
subprocess.call(["make", "debug"])
subprocess.call(["./main"])

# Find the error for runs 1-3 (as compared to run 4)
os.chdir("..")	
subprocess.call(["visit", "-cli","-s","getDiff.py"])


f = open('diff.txt','r')
diff = float(f.read())
f.close()

# Output the difference
print "Maximum difference between the calculated and analytic result: \n"
print str(diff), "\n"

text_file = open("diffusion_test_result.txt","w")
text_file.write("Maximum difference between the calculated and analytic result: \n") 
text_file.write(str(diff) + "\n") 
text_file.close()
