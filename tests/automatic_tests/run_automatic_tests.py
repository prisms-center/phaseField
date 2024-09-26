import subprocess
import shutil
import glob
import os
import datetime
import time
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed


def remove_file(file):
    """Remove a file if it exists

    Args:
            file (string): Path to file
    """
    if os.path.exists(file):
        os.remove(file)


def write_to_file(file, text):
    """Write text to file

    Args:
        file (string): Path to file
        text (string): Text to write in file
    """
    with open(file, "a") as text_file:
        text_file.write(text)


def compile_and_run_unit_tests():
    # Remove old files if they exist
    remove_file("main")
    remove_file("CMakeCache.txt")
    remove_file("unit_test_results.txt")

    # Compile and run
    subprocess.call(["cmake", "."])
    subprocess.call(["make"])
    subprocess.call(["mpirun", "-n", "1", "./main"])

    with open("unit_test_results.txt", "r") as result_file:
        return result_file.read().splitlines()


def compile_and_run_simulation(run_name):
    # Delete any pre-existing executables or results
    if os.path.exists(run_name):
        shutil.rmtree(run_name)

    # Remove old files
    remove_file("main")
    remove_file("CMakeCache.txt")
    for solution_file in glob.glob("solution-*.vtu"):
        os.remove(solution_file)

    # Compile application
    subprocess.call(["cmake", "."])
    subprocess.call(["make", "release"])

    # Print successful compilation to screen
    print("Compiling complete, running the regression test...")
    sys.stdout.flush()

    # Run application with timer
    start = time.time()
    subprocess.call(["mpirun", "-n", "1", "./main"])
    end = time.time()

    # Group the files
    subprocess.call(["mkdir", run_name])
    for output_files in glob.glob("solution-*.vtu"):
        shutil.move(output_files, run_name)
    if os.path.exists("integratedFields.txt"):
        shutil.move("integratedFields.txt", run_name)

    return end - start


def run_regression_test(applicationName, getNewGoldStandard, dir_path):

    testName = f"{'gold' if getNewGoldStandard else 'test'}_{applicationName}"

    if os.path.exists(testName):
        shutil.rmtree(testName)

    # Move to the application directory
    os.chdir(os.path.join("../../applications/", applicationName))

    # Run the simulation and move the results to the test directory
    test_time = compile_and_run_simulation(testName)

    shutil.move(testName, dir_path)
    os.chdir(dir_path)

    # Compare the result against the gold standard, if it exists
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

        rel_diff = abs(
            (float(gold_last_energy) - float(last_energy)) / float(gold_last_energy)
        )
        test_passed = rel_diff < 1.0e-9

    # Determine test result
    test_result = (
        "New Gold Standard"
        if test_passed and getNewGoldStandard
        else "Pass" if test_passed else "Fail"
    )

    # Print the results to the screen
    print(
        f"Regression Test: {applicationName}\n"
        f"Result: {test_result}\n"
        f"Time taken: {test_time}\n"
    )
    sys.stdout.flush()

    # Write the results to a file
    os.chdir(dir_path)
    write_to_file(
        "test_results.txt",
        f"Application: {applicationName}\n" f"Time: {test_time}\n \n",
    )

    return test_passed, test_time


def run_regression_tests_in_parallel(application_list, gold_standard_list, dir_path):

    regression_test_counter = 0
    regression_tests_passed = 0

    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(run_regression_test, application, gold_standard, dir_path)
            for application, gold_standard in zip(application_list, gold_standard_list)
        ]

        for future in as_completed(futures):
            test_passed, test_time = future.result()
            regression_tests_passed += int(test_passed)
            regression_test_counter += 1

    return regression_tests_passed, regression_test_counter


# Initialize
dir_path = os.path.dirname(os.path.realpath(__file__))
os.chdir(dir_path)
now = datetime.datetime.now()

os.chdir("../unit_tests/")
unit_test_results = compile_and_run_unit_tests()
unit_tests_passed = unit_test_results[0]
unit_test_counter = unit_test_results[1]

os.chdir(dir_path)
write_to_file(
    "test_results.txt",
    "--------------------------------------------------------- \n"
    "Unit test on " + now.strftime("%Y-%m-%d %H:%M") + "\n"
    "--------------------------------------------------------- \n"
    f"Unit Tests Passed: {unit_tests_passed}/{unit_test_counter}\n",
)


write_to_file(
    "test_results.txt",
    "--------------------------------------------------------- \n"
    "Regression test on " + now.strftime("%Y-%m-%d %H:%M") + "\n"
    "--------------------------------------------------------- \n",
)

# List of applications
applicationList = [
    "allenCahn",
    "cahnHilliard",
    "CHAC_anisotropyRegularized",
    "coupledCahnHilliardAllenCahn",
    "precipitateEvolution",
]
getNewGoldStandardList = [False, False, False, False, False]

regression_tests_passed, regression_test_counter = run_regression_tests_in_parallel(
    applicationList, getNewGoldStandardList, dir_path
)

# Output the overall test results
write_to_file(
    "test_results.txt",
    f"Tests Passed: {regression_tests_passed}/{regression_test_counter}\n"
    "--------------------------------------------------------- \n",
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
