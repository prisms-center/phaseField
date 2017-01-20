import subprocess
import shutil
import glob
import math
import os.path
import datetime

# ----------------------------------------------------------------------------------------
# Function that generates a header file for the simulation, compiles the PRISMS-PF code
# and runs the executable.
# ----------------------------------------------------------------------------------------
def run_simulation(run_name):

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

# Run simulation comparison simulation
run_simulation('run_001')

# Find the error for runs 1 compared to the reference solution
subprocess.call(["visit", "-cli","-s","getError.py"])

f = open('error.txt','r')
error = float(f.read())
f.close()

print "Regression test results: \n"
print "Difference between calculated result and reference result:", str(error), "\n"

text_file = open("regress_test_results.txt","a")
now = datetime.datetime.now()
text_file.write("Results of the 3D regression test (" + now.strftime("%Y-%m-%d %H:%M") + "): \n") 
text_file.write("Difference from reference: " + str(error) + "\n") 
if error < 1e-10:
	text_file.write("Result: Pass \n \n") 
else: 
	text_file.write("Result: Fail \n \n") 
text_file.close()

