import Rappture
#from Rappture.tools import executeCommand as RapptureExec
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

# ----------------------------------------------------------------------------------------
# Function that generates the PRISMS-PF input file snippet
# ----------------------------------------------------------------------------------------
def create_prismspf_input_file(entry_name, entry_value):

    if os.path.exists("parameters_rappture.in") is True:
        os.remove("parameters_rappture.in")

    shutil.copyfile("parameters.in", "parameters_rappture.in")

    f = open("parameters_rappture.in", 'a+')

    '''
    entry_counter = 0

    for name in entry_name:
        f.write("set " + name + " = " + entry_value[entry_counter] + "\n")
        entry_counter = entry_counter + 1
    '''

    for index, name in enumerate(entry_name):
        f.write("set " + entry_name[index] + " = " + entry_value[index] + "\n")

    f.close()


# ----------------------------------------------------------------------------------------
# Function that compiles the PRISMS-PF code and runs the executable.
# ----------------------------------------------------------------------------------------
def run_simulation(run_name, dir_path, num_time_steps, num_outputs):

    # Delete any pre-existing executables or results
    if os.path.exists(run_name) is True:
        shutil.rmtree(run_name)

    # Open file where output is redirected to
    if os.path.exists("output.txt") is True:
        os.remove("output.txt")
    f = open("output.txt", 'w+')

    for vtu_file in glob.glob("*.vtu"):
        os.remove(vtu_file)

    # Run the simulation
    #subprocess.call(["mpirun", "-n", "6", "main", "-i","parameters_rappture.in"], stdout=f)
    p = subprocess.Popen(["mpirun", "-n", "6", "main", "-i","parameters_rappture.in"], stdout=f)

    f = open('debug.txt','w')

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



    # Group the files
    subprocess.call(["mkdir", run_name])
    for output_files in glob.glob('*vtu'):
        shutil.move(output_files, run_name)
    if os.path.exists("integratedFields.txt") is True:
        shutil.move("integratedFields.txt", run_name)

# ----------------------------------------------------------------------------------------
# Function to extract a specific parameter from a PRISMS-PF input file
# ----------------------------------------------------------------------------------------
def parameter_extractor(file_name, entry_name):
    parameter_value = 0

    num_words_in_entry_name = len(entry_name.split())
    entry_name_no_whitespace = "".join(entry_name.split())

    f = open(file_name)
    for line in f:
        words_in_line = line.split()
        if len(words_in_line) > 0:
            if words_in_line[0] == "set":
                if (len(words_in_line) > num_words_in_entry_name+2):
                    prospective_match = ""
                    for word in range(1, 1+num_words_in_entry_name):
                        prospective_match += words_in_line[word]
                    if (prospective_match == entry_name_no_whitespace):
                        if len(words_in_line) >= 2+num_words_in_entry_name:
                            parameter_value = words_in_line[2+num_words_in_entry_name]

    f.close()
    return parameter_value

# ----------------------------------------------------------------------------------------
# Function to analyze the results using VisIt
# ----------------------------------------------------------------------------------------
def run_analysis(run_name, set_index, dir_path, num_time_steps):

    scratch_file = open("scratch.txt", 'w')
    scratch_file.write(str(set_index))
    scratch_file.write('\n')
    scratch_file.write(str(num_time_steps))
    scratch_file.close()

    subprocess.call(["visit", "-cli", "-nowin", "-s", "getPrecipDimAuto.py"])
    f = open("visit_analysis.txt", 'r')
    result = f.readlines()
    f.close()
    return result

# ----------------------------------------------------------------------------------------
# Function to calculate the minimum bounding rectangle from a set of points
# ----------------------------------------------------------------------------------------
def minimum_bounding_rectangle(points):
    """
    Find the smallest bounding rectangle for a set of points.
    Returns a set of points representing the corners of the bounding box.

    :param points: an nx2 matrix of coordinates
    :rval: an nx2 matrix of coordinates
    """
    from scipy.ndimage.interpolation import rotate
    pi2 = np.pi/2.

    # get the convex hull for the points
    hull_points = points[ConvexHull(points).vertices]

    # calculate edge angles
    edges = np.zeros((len(hull_points)-1, 2))
    edges = hull_points[1:] - hull_points[:-1]

    angles = np.zeros((len(edges)))
    angles = np.arctan2(edges[:, 1], edges[:, 0])

    angles = np.abs(np.mod(angles, pi2))
    angles = np.unique(angles)

    # find rotation matrices
    # XXX both work
    rotations = np.vstack([
        np.cos(angles),
        np.cos(angles-pi2),
        np.cos(angles+pi2),
        np.cos(angles)]).T
#     rotations = np.vstack([
#         np.cos(angles),
#         -np.sin(angles),
#         np.sin(angles),
#         np.cos(angles)]).T
    rotations = rotations.reshape((-1, 2, 2))

    # apply rotations to the hull
    rot_points = np.dot(rotations, hull_points.T)

    # find the bounding points
    min_x = np.nanmin(rot_points[:, 0], axis=1)
    max_x = np.nanmax(rot_points[:, 0], axis=1)
    min_y = np.nanmin(rot_points[:, 1], axis=1)
    max_y = np.nanmax(rot_points[:, 1], axis=1)

    # find the box with the best area
    areas = (max_x - min_x) * (max_y - min_y)
    best_idx = np.argmin(areas)

    # return the best box
    x1 = max_x[best_idx]
    x2 = min_x[best_idx]
    y1 = max_y[best_idx]
    y2 = min_y[best_idx]
    r = rotations[best_idx]

    rval = np.zeros((4, 2))
    rval[0] = np.dot([x1, y2], r)
    rval[1] = np.dot([x2, y2], r)
    rval[2] = np.dot([x2, y1], r)
    rval[3] = np.dot([x1, y1], r)

    return rval

# ----------------------------------------------------------------------
#  MAIN PROGRAM - generated by the Rappture Builder
# ----------------------------------------------------------------------

# open the XML file containing the run parameters
io = Rappture.library(sys.argv[1])

#########################################################
# Get input values from Rappture
#########################################################

# get input value for input.group(misfit_strains).number(misfit11)
interfacial_energy_11 = float(io.get('input.group(interfacial_energy).number(interfacial_energy_11).current'))
interfacial_energy_22 = float(io.get('input.group(interfacial_energy).number(interfacial_energy_22).current'))

# get input value for input.group(misfit_strains).number(misfit11)
misfit11 = float(io.get('input.group(misfit_strains).number(misfit11).current'))

# get input value for input.group(misfit_strains).number(misfit22)
misfit22 = float(io.get('input.group(misfit_strains).number(misfit22).current'))

# get input value for input.group(misfit_strains).number(misfit12)
misfit12 = float(io.get('input.group(misfit_strains).number(misfit12).current'))

# get input value for input.group(misfit_strains).number(misfit21)
misfit21 = float(io.get('input.group(misfit_strains).number(misfit21).current'))

# get input value for input.group(elastic_constants).group(ec_matrix).number(matrix_modulus)
matrix_modulus = float(io.get('input.group(ec_matrix).number(matrix_modulus).current'))

# get input value for input.group(elastic_constants).group(ec_matrix).number(matrix_poisson)
matrix_poisson = float(io.get('input.group(ec_matrix).number(matrix_poisson).current'))

# get input value for input.group(elastic_constants).group(ec_precip).number(precip_modulus)
precip_modulus = float(io.get('input.group(ec_precip).number(precip_modulus).current'))

# get input value for input.group(elastic_constants).group(ec_precip).number(precip_poisson)
precip_poisson = float(io.get('input.group(ec_precip).number(precip_poisson).current'))


#########################################################
#  Add your code here for the main body of your program
#########################################################

# Initialize
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

# spit out progress messages as you go along...
Rappture.Utils.progress(0, "Starting...")

#exitStatus,stdOutput,stdError = Rappture.tools.executeCommand('echo here')

misfit_string = '(('+str(misfit11)+','+str(misfit12)+',0),('+str(misfit21)+','+str(misfit22)+',0),(0,0,0)), tensor'

ec_matrix_string = '('+str(matrix_modulus/2.0)+','+str(matrix_poisson)+'), isotropic elastic constants'

ec_beta_string = '('+str(precip_modulus/2.0)+','+str(precip_poisson)+'), isotropic elastic constants'

interfacial_energy_string_11 = str(interfacial_energy_11)+', double'

interfacial_energy_string_22 = str(interfacial_energy_22)+', double'

create_prismspf_input_file(['Model constant sfts_const1', 'Model constant CIJ_Mg', 'Model constant CIJ_Beta', 'Model constant interfacial_energy_11', 'Model constant interfacial_energy_22'], (misfit_string, ec_matrix_string, ec_beta_string,interfacial_energy_string_11, interfacial_energy_string_22))


Rappture.Utils.progress(5, "Running the phase field simulation...")

num_time_steps = parameter_extractor("parameters_rappture.in", "Number of time steps")
num_outputs = parameter_extractor("parameters_rappture.in", "Number of outputs")

run_simulation("run_"+str(0), dir_path, float(num_time_steps), int(num_outputs))

Rappture.Utils.progress(90, "Simulation complete, beginning analysis...")



# Extract the points along the interface of the precipitate

scratch_file = open("scratch.txt", 'w')
scratch_file.write(str(0))
scratch_file.write('\n')
scratch_file.write(str(num_time_steps))
scratch_file.close()

#Rappture.Utils.progress(91, "Plotting the contours...")
subprocess.call(["visit", "-cli", "-nowin", "-s", "saveContour.py"])
#subprocess.call(["visit", "-cli", "-s", "saveContour.py"])

Rappture.Utils.progress(92, "Determining the precipitate dimensions...")

f = open('contour_ellipse.xyz','r')
point_list = []
for line in f:
    split_line = line.split()
    if len(split_line) > 4:
        point_list.append([float(split_line[1]), float(split_line[2])])

points = np.array(point_list)
rect = minimum_bounding_rectangle(points)

length_1 = sqrt((rect[0][0]-rect[1][0])**2 + (rect[0][1]-rect[1][1])**2)
length_2 = sqrt((rect[1][0]-rect[2][0])**2 + (rect[1][1]-rect[2][1])**2)

f = open("precipitate_lengths.txt",'w')
f.write(str(round(length_1,4))+'\n'+str(round(length_2,4)))
f.close()

#result = run_analysis("run_"+str(0), 0, dir_path, num_time_steps)

subprocess.call(["rm", "precipitate_plot_visit0000.png"])

Rappture.Utils.progress(93, "Plotting the precipitates...")
subprocess.call(["visit", "-cli", "-nowin", "-s", "plotPrecipitate.py"])
image = 'precipitate_plot_visit0000.png'

encoded_string = ""
with open(image, "rb") as image_file:
    encoded_string = base64.b64encode(image_file.read())
    imdata = encoded_string

aspect_ratio = max(length_1/length_2,length_2/length_1)

Rappture.Utils.progress(100, "Done")

#########################################################
# Save output values back to Rappture
#########################################################

# save output value for output.image(result_image)
# data should be base64-encoded image data
io.put('output.image(result_image).current', imdata)

# save output value for output.number(aspect_ratio)
io.put('output.number(aspect_ratio).current',aspect_ratio)

Rappture.result(io)
sys.exit()
