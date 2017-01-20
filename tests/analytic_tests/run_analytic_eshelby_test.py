import subprocess
import shutil
import glob
import math
import os.path
import os

# ----------------------------------------------------------------------------------------

# If files exist from previous tests, delete them
if os.path.exists("simulation_result_eshelby") == True:
	shutil.rmtree("simulation_result_eshelby")
if os.path.exists("analytical_result_eshelby") == True:
	shutil.rmtree("analytical_result_eshelby")

# Run the 1D diffusion simulation
os.chdir("eshelbyAnalyticComparison")
subprocess.call(["rm", "CMakeCache.txt"])
subprocess.call(["cmake", "."])
subprocess.call(["make", "release"])
subprocess.call(["./main"])

# Get the analytic result
os.chdir("../eshelbyAnalyticResult")
subprocess.call(["rm", "CMakeCache.txt"])
subprocess.call(["cmake", "."])
subprocess.call(["make", "release"])
subprocess.call(["./main"])

# Find the error for runs 1-3 (as compared to run 4)
os.chdir("..")	
subprocess.call(["visit", "-cli","-s","getDiff_eshelby.py"])


f = open('eshelby_difference.txt','r')
diff = float(f.read())
f.close()

# Output the difference
print "Maximum difference between the calculated and analytic result: \n"
print str(diff), "\n"

text_file = open("eshelby_test_result.txt","w")
text_file.write("Maximum difference between the calculated and analytic result: \n") 
text_file.write(str(diff) + "\n") 
text_file.close()

# Move the files
subprocess.call(["mkdir","simulation_result_eshelby"])
os.chdir("eshelbyAnalyticComparison")
subprocess.call(["cp *vtu ../simulation_result_eshelby/"],shell=True)

os.chdir("..")	
subprocess.call(["mkdir","analytical_result_eshelby"])
os.chdir("eshelbyAnalyticResult")
subprocess.call(["cp *vtu ../analytical_result_eshelby/"],shell=True)
os.chdir("..")	
