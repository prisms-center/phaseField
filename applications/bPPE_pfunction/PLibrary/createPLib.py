import os.path
import subprocess
import shutil
import glob
#import math
import os.path
import sys


# Get parameters that will be placed in the library (eventually through the MC API)
coeff_file = open('fa.txt','r')
fa_coeffs = coeff_file.read().splitlines()
coeff_file.close()

# Write PFunctions for the free energy and its first two derivatives 

f_writer_string = 'fw -n "pfunct_faV" -v "c" -d "concentration" --sym "'+fa_coeffs[0]+'*c^2 +'+fa_coeffs[1]+'*c + '+fa_coeffs[2]+'" --grad --hess -l "/Users/stephendewitt/Documents/workspace/PRISMS_workspace/PLibrary"'

subprocess.call([f_writer_string],shell=True)

l_writer_string = 'lw -d /Users/stephendewitt/Documents/workspace/PRISMS_workspace/PLibrary -v "std::vector<double>" -l /Users/stephendewitt/Documents/workspace/PRISMS_workspace/PLibrary -c'

subprocess.call([l_writer_string],shell=True)


