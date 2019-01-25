import Rappture
from Rappture.tools import executeCommand as RapptureExec
import sys
from math import *

import subprocess
import shutil
import glob
import os
import datetime
import time

import base64

import numpy as np
from scipy.spatial import ConvexHull

if __name__ == "__main__":
    
    f = open('debug.txt','w')
    p = subprocess.Popen(["mpirun", "-n", "6", "main", "-i","parameters_rappture.in"], stdout=f)

    fts = open("run_info.txt")
    num_time_steps = int(fts.readline().rstrip())
    num_outputs = int(fts.readline())
    fts.close()
    
    # Print out progress during the simulation
    solution_file_names = []
    f.write(str(num_time_steps))
    f.close()

    for findex in range(1, num_outputs + 1):
        time_step_number = int(num_time_steps / num_outputs * findex)
        time_step_number_string = str(time_step_number).rjust(len(str(int(num_time_steps))), '0')
        name = 'solution-' + time_step_number_string + '.vtu'
        solution_file_names.append(name)
        f = open('debug.txt', 'a')
        f.write(name + '\n')
        f.close()

    output_index = 0
    converged = False
    while True:
        time.sleep(0.5)
        f = open('debug.txt', 'a')
        f.write(solution_file_names[output_index] + '\n')
        f.close()
        
        if os.path.exists(solution_file_names[output_index]):
            progress = int((85.0 - 5.0) / num_outputs * (output_index + 1) + 5)
            Rappture.Utils.progress(progress, "Running the phase field simulation...")

            # Check to see if the free energy has converged
            energy_file = open('integratedFields.txt', 'r')
            energy_file_contents = energy_file.readlines()
            current_energy = energy_file_contents[output_index + 1].split()[2]
            prior_energy = energy_file_contents[output_index].split()[2]
            f = open('debug.txt', 'a')
            f.write(str(current_energy)+ ' ' + str(prior_energy) + '\n')
            f.close()
            if (abs(float(current_energy) - float(prior_energy)) < 1e-9):
	        f = open('debug.txt', 'a')
	        f.write('killed! \n')
	        f.close()

	        converged = True
	        p.kill()

	    output_index = output_index + 1

	if (output_index >= num_outputs) or converged:
	    break

    f = open('debug.txt', 'a')
    f.write('Out of check loop \n')
    f.close()

    if converged:
	for findex in range(output_index, num_outputs):
	    shutil.copyfile(solution_file_names[output_index - 1], solution_file_names[findex])

    p.wait()
