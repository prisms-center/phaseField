import os.path
import subprocess
import shutil
import glob
#import math
import os.path
import sys

# -----------------------------------------------------------------------------------------
# Function that creates a PFunction using the IntegrationToolsWriter command line interface
# -----------------------------------------------------------------------------------------
# Inputs:
# f_name 	 = The name the function is referred to in the code that calls the PLibrary
# var		 = The variable(s) used in the function express, should match the variable names in the code that calls the PLibrary
# expression = The expression for the desired function in terms of "var"
# need_grad	 = Flag set to true if the gradient of the expression is needed, false if not
# need_hess	 = Flag set to true if the hessian of the expression is needed, false if not
# dir		 = Directory to write the PFunction to

def write_pfunction(f_name, var, var_name, expression, need_grad, need_hess, dir):
	f_writer_string = 'fw -n ' + f_name + ' -v ' + var + ' -d ' + var_name + ' --sym "'+ expression +'"'
	if need_grad == True:
		f_writer_string += " --grad "
	if need_hess == True:
		f_writer_string += " --hess "
		
	f_writer_string += ' -l '+dir
	
	subprocess.call([f_writer_string],shell=True)

# -----------------------------------------------------------------------------------------
# Function that creates a PLibrary using the IntegrationTools command line interface
# -----------------------------------------------------------------------------------------
# Inputs:
# data_type	 	= Data type for the container holding the PLibrary
# pfunction_dir = Directory containing the PFunctions to be assembled into the PLibrary
# plibrary_dir 	= Directory where the PLibrary should be written

def write_plibrary(data_type, pfunction_dir, plibrary_dir):
	l_writer_string = 'lw -d ' + pfunction_dir + ' -v "'+ data_type + '" -l '+ plibrary_dir +' -c --include "<deal.II/base/vectorization.h>"'
	
	print l_writer_string
	
	subprocess.call([l_writer_string],shell=True)

# -----------------------------------------------------------------------------------------
# Main script to generate a PLibrary
# -----------------------------------------------------------------------------------------

# Get parameters that will be placed in the library (eventually through the MC API, currently in text files)
coeff_file = open('fa.txt','r')
fa_coeffs = coeff_file.read().splitlines()
coeff_file.close()

coeff_file = open('fb.txt','r')
fb_coeffs = coeff_file.read().splitlines()
coeff_file.close()

CH_mobility_file = open('Mc.txt','r')
Mc = CH_mobility_file.read()
CH_mobility_file.close()

AC_mobility_file = open('Mn.txt','r')
Mn = AC_mobility_file.read().splitlines()
AC_mobility_file.close()

# Get the current directory, which is where the PFunctions and PLibrary will be created
dir = os.getcwd()

# Write PFunctions for the free energies and their first two derivatives 
write_pfunction("pfunct_faV", "c", "concentration", fa_coeffs[0]+'*c^4 +'+fa_coeffs[1]+'*c^3 + '+fa_coeffs[2]+'*c^2 +'+fa_coeffs[3]+'*c +'+fa_coeffs[4], True, True, dir)
write_pfunction("pfunct_fbV", "c", "concentration", fb_coeffs[0]+'*c^2 +'+fb_coeffs[1]+'*c +'+fb_coeffs[2], True, True, dir)

# Write a PFunction for the Cahn-Hilliard mobility
write_pfunction("pfunct_McV", "c", "concentration", Mc, False, False, dir)

# Write PFunctions for the Allen-Cahn mobilities
write_pfunction("pfunct_Mn1V", "n1", "concentration", Mn[0], False, False, dir)
write_pfunction("pfunct_Mn2V", "n2", "concentration", Mn[1], False, False, dir)
write_pfunction("pfunct_Mn3V", "n3", "concentration", Mn[2], False, False, dir)

# Write the PLibrary
write_plibrary("dealii::VectorizedArray<double>", dir, dir)



