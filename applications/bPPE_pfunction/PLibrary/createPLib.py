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

CH_mobility_file = open('Mc.txt','r')
Mc = CH_mobility_file.read()
CH_mobility_file.close()

# Write a PFunction for the free energy and its first two derivatives 
dir = os.getcwd()

f_writer_string = 'fw -n "pfunct_faV" -v "c" -d "concentration" --sym "'+fa_coeffs[0]+'*c^2 +'+fa_coeffs[1]+'*c + '+fa_coeffs[2]+'" --grad --hess -l '+dir

subprocess.call([f_writer_string],shell=True)

# Write a PFunction for the Cahn-Hilliard Mobility
f_writer_string = 'fw -n "pfunct_McV"  -v "c" -d "concentration" --sym "'+Mc+'" -l '+dir

subprocess.call([f_writer_string],shell=True)

# Write the PLibrary
l_writer_string = 'lw -d '+dir+' -v "dealii::VectorizedArray<double>" -l '+dir+' -c --include "<deal.II/base/vectorization.h>"'

subprocess.call([l_writer_string],shell=True)


