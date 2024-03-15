import subprocess
import shutil
import glob
import math
import os
import datetime
import time
import sys

# ----------------------------------------------------------------------------------------
# Function that compiles and runs the unit tests.
# ----------------------------------------------------------------------------------------
def run_unit_tests():
	# Open file where output is redirected to
	if os.path.exists("output.txt") == True:
		os.remove("output.txt")
	f = open("output.txt",'w+')

	# Remove old files
	if os.path.exists("main") == True:
		os.remove("main")
	if os.path.exists("CMakeCache.txt") == True:
		os.remove("CMakeCache.txt")
	if os.path.exists("unit_test_results.txt") == True:
		os.remove("unit_test_results.txt")

	# Compile and run
	subprocess.call(["cmake", "."],stdout=f,stderr=f)
	subprocess.call(["make", "release"],stdout=f)
	subprocess.call(["make"],stdout=f)
	subprocess.call(["mpirun", "-n", "2", "./main"],stdout=f)
	f.close()

	result_file = open("unit_test_results.txt","r")
	test_results = result_file.read().splitlines()
	result_file.close()

	return test_results


# ----------------------------------------------------------------------------------------
# Function that compiles the PRISMS-PF code and runs the executable.
# ----------------------------------------------------------------------------------------
def run_simulation(run_name,dir_path):
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
	subprocess.call(["cmake", "."]) # print to the screen to prevent a timeout on Travis
	subprocess.call(["make", "release"],stdout=f)
	print("Compiling complete, running the regression test...")
	sys.stdout.flush()
	start = time.time()
	subprocess.call(["mpirun", "-n", "1", "./main"],stdout=f)
	end = time.time()
	f.close()

	# Group the files
	subprocess.call(["mkdir",run_name])
	for output_files in glob.glob('*vtu'):
		shutil.move(output_files,run_name)
	if os.path.exists("integratedFields.txt") == True:
		shutil.move("integratedFields.txt",run_name)

	test_time = end-start
	return test_time

# ----------------------------------------------------------------------------------------
# Function that runs a regression test
# ----------------------------------------------------------------------------------------
def run_regression_test(applicationName,getNewGoldStandard,dir_path):

	if (getNewGoldStandard == False):
		testName = "test_"+applicationName

	else:
		testName = "gold_"+applicationName

	if os.path.exists(testName) == True:
		shutil.rmtree(testName)

	# Move to the application directory
	r_test_dir = dir_path
	os.chdir("../../applications/"+applicationName)

	# Run the simulation and move the results to the test directory
	test_time = run_simulation(testName,dir_path)

	shutil.move(testName,r_test_dir)

	# Compare the result against the gold standard, if it exists
	os.chdir(r_test_dir)

	if (getNewGoldStandard == False):
		# Read the gold standard free energies
		os.chdir("gold_"+applicationName)
		gold_standard_file = open("integratedFields.txt","r")
		gold_energy = gold_standard_file.readlines()
		gold_standard_file.close()

		last_energy_index = len(gold_energy)-1
		split_last_line = gold_energy[-1].split()
		for index, entry in enumerate(split_last_line):
			if entry == "f_tot":
				gold_last_energy = split_last_line[index+1]

		# Read the test free energies
		os.chdir("../"+testName)
		test_file = open("integratedFields.txt","r")
		test_energy = test_file.readlines()
		test_file.close()

		last_energy_index = len(test_energy)-1
		split_last_line = test_energy[-1].split()
		for index, entry in enumerate(split_last_line):
			if entry == "f_tot":
				last_energy = split_last_line[index+1]

		rel_diff = (float(gold_last_energy)-float(last_energy))/float(gold_last_energy)
		rel_diff = abs(rel_diff)

		if (rel_diff < 1.0e-9):
			test_passed = True
		else:
			test_passed = False

	else:
		test_passed = True

	# Print the results to the screen
	print("Regression Test: ", applicationName)

	if test_passed:
		if getNewGoldStandard == False:
			print("Result: Pass")
		else:
			print("Result: New Gold Standard")
	else:
		print("Result: Fail")

	print("Time taken:", test_time)

	sys.stdout.flush()

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
# Test Script
# ----------------------------------------------------------------------------------------

# Initialize
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

text_file = open("test_results.txt","a")
now = datetime.datetime.now()

# ------------------------------------
# Start the unit tests
# ------------------------------------


os.chdir("../unit_tests/")
unit_test_results = run_unit_tests()
unit_tests_passed = unit_test_results[0]
unit_test_counter = unit_test_results[1]

print()
print("Unit Tests Passed: "+str(unit_tests_passed)+"/"+str(unit_test_counter)+"\n")

sys.stdout.flush()

text_file.write("--------------------------------------------------------- \n")
text_file.write("Unit test on " + now.strftime("%Y-%m-%d %H:%M") + "\n")
text_file.write("--------------------------------------------------------- \n")
text_file.write("Unit Tests Passed: "+str(unit_tests_passed)+"/"+str(unit_test_counter)+"\n")

os.chdir(dir_path)

# ------------------------------------
# Start the regression tests
# ------------------------------------

regression_test_counter = 0
regression_tests_passed = 0

text_file.write("--------------------------------------------------------- \n")
text_file.write("Regression test on " + now.strftime("%Y-%m-%d %H:%M") + "\n")
text_file.write("--------------------------------------------------------- \n")
text_file.close()

#applicationList = ["allenCahn","cahnHilliard","cahnHilliard","CHAC_anisotropy","CHAC_anisotropyRegularized","coupledCahnHilliardAllenCahn","mechanics","precipitateEvolution"]
#getNewGoldStandardList = [False, False, False, False, False, False, False, False]

# Shorter list of applications so that it completes on Travis
applicationList = ["allenCahn","cahnHilliard","CHAC_anisotropyRegularized","coupledCahnHilliardAllenCahn","precipitateEvolution"]
getNewGoldStandardList = [False, False, False, False, False]
#applicationList = ["allenCahn"]
#getNewGoldStandardList = [False]



for applicationName in applicationList:
	test_result = run_regression_test(applicationName,getNewGoldStandardList[regression_test_counter],dir_path)

	regression_test_counter += 1
	regression_tests_passed += int(test_result[0])

print()
print("Regression Tests Passed: "+str(regression_tests_passed)+"/"+str(regression_test_counter)+"\n")


# Output the overall test results
text_file = open("test_results.txt","a")
text_file.write("Tests Passed: "+str(regression_tests_passed)+"/"+str(regression_test_counter)+"\n")
text_file.write("--------------------------------------------------------- \n")
text_file.close()

# Set exit code (passed to Travis CI)
if ((regression_tests_passed < regression_test_counter) or (unit_tests_passed < unit_test_counter)):
	sys.exit(1)
else:
	sys.exit(0)
