import subprocess
import shutil
import glob
import math
import os.path
import os

# ----------------------------------------------------------------------------------------

# If files exist from previous tests, delete them
if os.path.exists("simulation_result_diffusion") == True:
	shutil.rmtree("simulation_result_diffusion")
if os.path.exists("analytical_result_diffusion") == True:
	shutil.rmtree("analytical_result_diffusion")

# Run the 1D diffusion simulation

os.chdir("diffusionAnalyticComparison")
subprocess.call(["rm", "CMakeCache.txt"])
subprocess.call(["cmake", "."])
subprocess.call(["make", "debug"])
subprocess.call(["./main"])

# Get the analytic result
os.chdir("../diffusionAnalyticResult")
subprocess.call(["rm", "CMakeCache.txt"])
subprocess.call(["cmake", "."])
subprocess.call(["make", "debug"])
subprocess.call(["./main"])

# Find the error 
os.chdir("..")	
subprocess.call(["visit", "-cli","-s","getDiff_diffusion.py"])

f = open('diffusion_difference.txt','r')
diff = float(f.read())
f.close()

# Output the difference
print "Maximum difference between the calculated and analytic result: \n"
print str(diff), "\n"

text_file = open("diffusion_test_result.txt","w")
text_file.write("Maximum difference between the calculated and analytic result: \n") 
text_file.write(str(diff) + "\n") 
text_file.close()

# Move the files
subprocess.call(["mkdir","simulation_result_diffusion"])
os.chdir("diffusionAnalyticComparison")
subprocess.call(["cp *vtu ../simulation_result_diffusion/"],shell=True)

os.chdir("..")	
subprocess.call(["mkdir","analytical_result_diffusion"])
os.chdir("diffusionAnalyticResult")
subprocess.call(["cp *vtu ../analytical_result_diffusion/"],shell=True)
os.chdir("..")	

