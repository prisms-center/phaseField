import subprocess
import shutil
import glob
import math
import os
import datetime
import time
import sys


# ----------------------------------------------------------------------------------------
# Function that generates the PRISMS-PF input file snippet
# ----------------------------------------------------------------------------------------
def create_prismspf_input_file(entry_name, entry_value):

    # if os.path.exists("additional_parameters.in") is True:
    #     os.remove("additional_parameters.in")
    # f = open("additional_parameters.in", 'w+')
    #
    # entry_counter = 0
    #
    # for name in entry_name:
    #     f.write("set " + name + " = " + entry_value[entry_counter] + ", double\n")
    #     entry_counter += 1
    #
    # f.close()

    if os.path.exists("parameters.in") is True:
        os.remove("parameters.in")

    shutil.copyfile("parameters_base.in", "parameters.in")

    f = open("parameters.in", 'a+')

    entry_counter = 0

    for name in entry_name:
        f.write("set " + name + " = " + entry_value[entry_counter] + ", double\n")
        entry_counter += 1

    f.close()


# ----------------------------------------------------------------------------------------
# Function that compiles the PRISMS-PF code and runs the executable.
# ----------------------------------------------------------------------------------------
def run_simulation(run_name, dir_path):

    # Delete any pre-existing executables or results
    if os.path.exists(run_name) is True:
        shutil.rmtree(run_name)

    # Open file where output is redirected to
    if os.path.exists("output.txt") is True:
        os.remove("output.txt")
    f = open("output.txt", 'w+')

    subprocess.call(["rm", "*vtu"], stdout=f, stderr=f)

    # Run the simulation
    subprocess.call(["mpirun", "-n", "2", "main"], stdout=f)

    # Group the files
    subprocess.call(["mkdir", run_name])
    for output_files in glob.glob('*vtu'):
        shutil.move(output_files, run_name)
    if os.path.exists("integratedFields.txt") is True:
        shutil.move("integratedFields.txt", run_name)


def run_analysis(run_name, dir_path):
    subprocess.call(["visit", "-cli", "-nowin", "-s", "getPrecipDimAuto.py"])
    f = open("visit_analysis.txt", 'r')
    result = f.readlines()
    f.close()
    return result


def write_dakota_file(run_num, result):
    f = open("results.out."+str(run_num), 'w+')
    for line in result:
        f.write(line.strip() + " f" + "\n")
    f.close()


def read_dakota_files():
    # Determine how many files need to be read
    max_index = 1

    while os.path.exists("params.in."+str(max_index)):
        max_index += 1

    max_index -= 1

    print "max_index", max_index

    # Load the variable entry names and values from the files
    parameter_set_list = []
    for index in range(1, max_index+1):
        f = open("params.in."+str(index), 'r')
        file_contents = f.readlines()
        f.close()
        temp = file_contents[0].replace("variables", "")
        num_parameters = int(temp.strip())

        print "num parameters: ", num_parameters

        entry_name = []
        entry_value = []
        for line in range(1, num_parameters+1):
            temp = file_contents[line].strip()
            split_strings = temp.split()

            entry_value.append(split_strings[0])
            entry_name_temp = ""
            for item in split_strings[1:len(split_strings)]:
                entry_name_temp += " " + item
                print entry_name_temp
            entry_name.append(entry_name_temp)

            print "parameter set", entry_name, " ", entry_value

        parameter_set_list.append((entry_name, entry_value))

    return parameter_set_list

# ----------------------------------------------------------------------------------------
# Dakota Interface Script
# ----------------------------------------------------------------------------------------
# Initialize
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

# Read input from Dakota
parameter_set_list = read_dakota_files()

# ------------------------------------
# Run the PRISMS-PF App
# ------------------------------------
print parameter_set_list

set_index = 1
for parameter_set in parameter_set_list:
    result = ["0.0", "0.0", "0.0"]

    write_dakota_file(set_index, result)
    set_index += 1
