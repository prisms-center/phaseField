import subprocess
import shutil
import glob
import math
import os
import datetime
import time
import sys

# ----------------------------------------------------------------------------------------
# Function tha compiles the PRISMS-PF code and runs the executable.
# ----------------------------------------------------------------------------------------
def run_simulation(run_name):
	# Delete any pre-existing executables or results
	if os.path.exists(run_name) == True:
		shutil.rmtree(run_name)
	
	# Open file where output is redirected to
	if os.path.exists("output.txt") == True:
		os.remove("output.txt")
	f = open("output.txt",'w+')	
	
	# Remove old files
	if os.path.exists("main") == True:
		os.remove("main")
	if os.path.exists("CMakeCache.txt") == True:
		os.remove("CMakeCache.txt")
		
	subprocess.call(["rm", "*vtu"],stdout=f,stderr=f)
	
	# Compile and run
	subprocess.call(["cmake", "."],stdout=f,stderr=f)
	subprocess.call(["make", "release"],stdout=f)
	start = time.time()
	subprocess.call(["mpirun", "-n", "4", "main"],stdout=f)
	end = time.time()
	f.close()

	# Group the files
	subprocess.call(["mkdir",run_name])
	for output_files in glob.glob('*vtu'):
		shutil.move(output_files,run_name)
	if os.path.exists("freeEnergy.txt") == True:
		shutil.move("freeEnergy.txt",run_name)
	
	test_time = end-start
	return test_time

# ----------------------------------------------------------------------------------------
# Function that runs a regression test
# ----------------------------------------------------------------------------------------
def run_regression_test(applicationName,getNewGoldStandard):	

	if (getNewGoldStandard == False):
		testName = "test_"+applicationName
	
	else:
		testName = "gold_"+applicationName

	if os.path.exists(testName) == True:
		shutil.rmtree(testName)

	# Move to the application directory
	r_test_dir = os.getcwd()
	os.chdir("../../applications/"+applicationName)

	# Run the simulation and move the results to the test directory
	test_time = run_simulation(testName)

	shutil.move(testName,r_test_dir)

	# Compare the result against the gold standard, if it exists
	os.chdir(r_test_dir)

	if (getNewGoldStandard == False):
		# Read the gold standard free energies
		os.chdir("gold_"+applicationName)
		gold_standard_file = open("freeEnergy.txt","r")
		gold_energy = gold_standard_file.readlines()
		gold_standard_file.close()
	
		# Read the test free energies
		os.chdir("../"+testName)
		test_file = open("freeEnergy.txt","r")
		test_energy = test_file.readlines()
		test_file.close()
	
		rel_diff = (float(gold_energy[9])-float(test_energy[9]))/float(gold_energy[9])
		rel_diff = abs(rel_diff)
	
		if (rel_diff < 1.0e-10):
			test_passed = True
		else:
			test_passed = False
		
	else:
		test_passed = True
		
	# Print the results to the screen
	print "Test: ", applicationName
	print "Passed: ", test_passed
	print "Time taken:", test_time
	
	# Write the results to a file
	os.chdir(r_test_dir)
	text_file = open("test_results.txt","a")
	now = datetime.datetime.now()
	text_file.write("Application: " + applicationName +" \n") 
	if test_passed:
		if getNewGoldStandard == False:
			text_file.write("Result: Pass \n") 
		else:
			text_file.write("Result: New Gold Standard \n") 
	else: 
		text_file.write("Result: Fail \n") 
	text_file.write("Time: "+str(test_time)+" \n \n") 
	text_file.close()
	
	return (test_passed,test_time)
# ----------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------


# Initialize
test_counter = 0
tests_passed = 0

text_file = open("test_results.txt","a")
now = datetime.datetime.now()
text_file.write("--------------------------------------------------------- \n")
text_file.write("Regression test on " + now.strftime("%Y-%m-%d %H:%M") + "\n") 
text_file.write("--------------------------------------------------------- \n")
text_file.close()

# Start the tests
# applicationName = "cahnHilliard"
# getNewGoldStandard = False
# 
# test_result = run_regression_test(applicationName,getNewGoldStandard)	
# 
# test_counter += 1
# tests_passed += int(test_result[0])

applicationName = "allenCahn"
getNewGoldStandard = False

test_result = run_regression_test(applicationName,getNewGoldStandard)	

test_counter += 1
tests_passed += int(test_result[0])

# Output the overall test results
text_file = open("test_results.txt","a")
text_file.write("Tests Passed: "+str(tests_passed)+"/"+str(test_counter)+"\n") 
text_file.write("--------------------------------------------------------- \n")
text_file.close()

if tests_passed < test_counter:
	sys.exit(1)
else:
	sys.exit(0)
	


