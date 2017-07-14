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
    subprocess.call(["mpirun", "-n", "4", "main"], stdout=f)

    # Group the files
    subprocess.call(["mkdir", run_name])
    for output_files in glob.glob('*vtu'):
        shutil.move(output_files, run_name)
    if os.path.exists("integratedFields.txt") is True:
        shutil.move("integratedFields.txt", run_name)


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

        entry_name = []
        entry_value = []
        for line in range(1, num_parameters+1):
            temp = file_contents[line].strip()
            split_strings = temp.split()

            entry_value.append(split_strings[0])
            entry_name_temp = ""
            for item in split_strings[1:len(split_strings)]:
                entry_name_temp += " " + item

            entry_name.append(entry_name_temp)

        parameter_set_list.append((entry_name, entry_value))

    return parameter_set_list


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
# Dakota Interface Script
# ----------------------------------------------------------------------------------------
# Initialize
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

subprocess.call(["rm", "run_*"])

# Read input from Dakota
run_num = 1
parameter_set_list = read_dakota_files()

# Compile
subprocess.call(["cmake", "."])
subprocess.call(["make", "release", "-j9"])
print "Compiling complete, running simulations now..."

# ------------------------------------
# Run the PRISMS-PF App
# ------------------------------------

set_index = 1
for parameter_set in parameter_set_list:

    entry_name = parameter_set[0]
    entry_value = parameter_set[1]

    create_prismspf_input_file(entry_name, entry_value)
    num_time_steps = parameter_extractor("parameters.in", "Number of time steps")
    run_simulation("run_"+str(set_index), dir_path)
    result = run_analysis("run_"+str(set_index), set_index, dir_path, num_time_steps)
    write_dakota_file(set_index, result)
    set_index += 1
