import subprocess
import shutil
import glob
import os
import datetime
import time
import sys
from concurrent.futures import ProcessPoolExecutor


def remove_file(file):
    """Remove a file if it exists

    Args:
            file (string): path of the file
    """
    if os.path.exists(file):
        os.remove(file)


def run_unit_tests():
    # Remove old files if they exist
    remove_file("output.txt")
    remove_file("main")
    remove_file("CMakeCache.txt")
    remove_file("unit_test_results.txt")

    # Open file where output is redirected to
    f = open("output.txt", "w+")

    # Compile and run
    subprocess.call(["cmake", "."], stdout=f, stderr=f)
    subprocess.call(["make"], stdout=f)
    subprocess.call(["mpirun", "-n", "1", "./main"], stdout=f)
    f.close()

    result_file = open("unit_test_results.txt", "r")
    test_results = result_file.read().splitlines()
    result_file.close()

    return test_results


def run_simulation(run_name):
    # Delete any pre-existing executables or results
    if os.path.exists(run_name):
        shutil.rmtree(run_name)

    # Open file where output is redirected to
    remove_file("output.txt")
    f = open("output.txt", "w+")

    # Remove old files
    if os.path.exists("main"):
        os.remove("main")
    if os.path.exists("CMakeCache.txt"):
        os.remove("CMakeCache.txt")

    subprocess.call(["rm", "*vtu"], stdout=f, stderr=f)

    # Compile and run
    subprocess.call(["cmake", "."])
    subprocess.call(["make", "release"], stdout=f)
    print("Compiling complete, running the regression test...")
    sys.stdout.flush()
    start = time.time()
    subprocess.call(["mpirun", "-n", "1", "./main"], stdout=f)
    end = time.time()
    f.close()

    # Group the files
    subprocess.call(["mkdir", run_name])
    for output_files in glob.glob("*vtu"):
        shutil.move(output_files, run_name)
    if os.path.exists("integratedFields.txt"):
        shutil.move("integratedFields.txt", run_name)

    test_time = end - start
    return test_time


def run_regression_test(applicationName, getNewGoldStandard, dir_path):

    if not getNewGoldStandard:
        testName = "test_" + applicationName
    else:
        testName = "gold_" + applicationName

    if os.path.exists(testName):
        shutil.rmtree(testName)

    # Move to the application directory
    r_test_dir = dir_path
    os.chdir("../../applications/" + applicationName)

    # Run the simulation and move the results to the test directory
    test_time = run_simulation(testName)

    shutil.move(testName, r_test_dir)

    # Compare the result against the gold standard, if it exists
    os.chdir(r_test_dir)

    if getNewGoldStandard:
        test_passed = True
    else:
        # Read the gold standard free energies
        os.chdir("gold_" + applicationName)
        gold_standard_file = open("integratedFields.txt", "r")
        gold_energy = gold_standard_file.readlines()
        gold_standard_file.close()

        split_last_line = gold_energy[-1].split()
        for index, entry in enumerate(split_last_line):
            if entry == "f_tot":
                gold_last_energy = split_last_line[index + 1]

        # Read the test free energies
        os.chdir("../" + testName)
        test_file = open("integratedFields.txt", "r")
        test_energy = test_file.readlines()
        test_file.close()

        split_last_line = test_energy[-1].split()
        for index, entry in enumerate(split_last_line):
            if entry == "f_tot":
                last_energy = split_last_line[index + 1]

        rel_diff = (float(gold_last_energy) - float(last_energy)) / float(
            gold_last_energy
        )
        rel_diff = abs(rel_diff)

        if rel_diff < 1.0e-9:
            test_passed = True
        else:
            test_passed = False

    # Print the results to the screen
    print(f"Regression Test: {applicationName}")

    if test_passed:
        if getNewGoldStandard:
            print("Result: New Gold Standard")
        else:
            print("Result: Pass")
    else:
        print("Result: Fail")

    print(f"Time taken: {test_time}")

    sys.stdout.flush()

    # Write the results to a file
    os.chdir(r_test_dir)
    with open("test_results.txt", "a") as text_file:
        text_file.write(f"Application: {applicationName}\n")
        if test_passed:
            if getNewGoldStandard:
                text_file.write("Result: New Gold Standard \n")
            else:
                text_file.write("Result: Pass \n")
        else:
            text_file.write("Result: Fail \n")
        text_file.write(f"Time: {test_time}\n \n")

    return (test_passed, test_time)


# Initialize
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)

now = datetime.datetime.now()


os.chdir("../unit_tests/")
unit_test_results = run_unit_tests()
unit_tests_passed = unit_test_results[0]
unit_test_counter = unit_test_results[1]


sys.stdout.flush()
with open("test_results.txt", "a") as text_file:
    text_file.write(
        "--------------------------------------------------------- \n"
        "Unit test on " + now.strftime("%Y-%m-%d %H:%M") + "\n"
        "--------------------------------------------------------- \n"
        f"Unit Tests Passed: {unit_tests_passed}/{unit_test_counter}\n"
    )

os.chdir(dir_path)


with open("test_results.txt", "a") as text_file:
    text_file.write(
        "--------------------------------------------------------- \n"
        "Regression test on " + now.strftime("%Y-%m-%d %H:%M") + "\n"
        "--------------------------------------------------------- \n"
    )

# Shorter list of applications so that it completes on Travis
applicationList = [
    "allenCahn",
    "cahnHilliard",
    "CHAC_anisotropyRegularized",
    "coupledCahnHilliardAllenCahn",
    "precipitateEvolution",
]
getNewGoldStandardList = [False, False, False, False, False]

regression_test_counter = 0
regression_tests_passed = 0

for applicationName in applicationList:
    test_result = run_regression_test(
        applicationName, getNewGoldStandardList[regression_test_counter], dir_path
    )

    regression_test_counter += 1
    regression_tests_passed += int(test_result[0])


# Output the overall test results
with open("test_results.txt", "a") as text_file:
    text_file.write(
        f"Tests Passed: {regression_tests_passed}/{regression_test_counter}\n"
        "--------------------------------------------------------- \n"
    )

# Print overall results
print(f"Unit Tests Passed: {unit_tests_passed}/{unit_test_counter}\n")
print(f"Regression Tests Passed: {regression_tests_passed}/{regression_test_counter}\n")

# Set exit code
if (regression_tests_passed < regression_test_counter) or (
    unit_tests_passed < unit_test_counter
):
    sys.exit(1)
else:
    sys.exit(0)
