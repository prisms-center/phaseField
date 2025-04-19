import subprocess
import shutil
import glob
import os
import datetime
import time
import sys
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import re


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


def grab_cpu_information():
    # Get the CPU information
    cpu_info_raw = subprocess.check_output("lscpu", shell=True).decode()

    # Split into lines
    cpu_info_lines = cpu_info_raw.splitlines()

    # Create dict to store parsed info
    cpu_info = {}

    for line in cpu_info_lines:
        if ":" in line:
            key, value = line.split(":", 1)
            cpu_info[key.strip()] = value.strip()

    # Grab relevant fields
    architecture = cpu_info.get("Architecture")
    cpu_model = cpu_info.get("Model name")
    cpu_cores = cpu_info.get("CPU(s)")
    cpu_max_freq = cpu_info.get("CPU max MHz")
    cpu_min_freq = cpu_info.get("CPU min MHz")
    hypervisor = cpu_info.get("Hypervisor vendor")

    return architecture, cpu_model, cpu_cores, cpu_max_freq, cpu_min_freq, hypervisor


def compile_and_run_simulation(application_path, n_threads=1):
    # Navigate to test application directory
    os.chdir(application_path)

    # Remove old files
    remove_file("main")
    remove_file("main-debug")
    remove_file("CMakeCache.txt")
    for solution_file in glob.glob("solution-*.vtu"):
        os.remove(solution_file)
    for solution_file in glob.glob("solution-*.pvtu"):
        os.remove(solution_file)

    # Compile application
    try:
        subprocess.run(
            ["cmake", ".", "-G", "Ninja"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=True,
        )
        subprocess.run(
            ["ninja", "-j", f"{n_threads}"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=True,
        )
        print("Compiling complete.")
    except subprocess.CalledProcessError as e:
        # Print the error stream if there's a compilation error
        print(
            f"Compilation failed with error:\n{e.stderr.decode('utf-8')}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Print successful compilation to screen
    print(f"Running the regression test for {application_path}...")
    sys.stdout.flush()

    # Run application with timer
    start = time.time()
    with open("output.txt", "w") as outfile:
        try:
            subprocess.run(
                ["mpirun", "-n", f"{n_threads}", "./main"],
                stdout=outfile,
                stderr=subprocess.PIPE,
                check=True,
            )
        except subprocess.CalledProcessError as e:
            # Print the error stream if the application fails
            print(
                f"Application execution failed with error:\n{e.stderr.decode('utf-8')}",
                file=sys.stderr,
            )
            sys.exit(1)

    end = time.time()
    print(f"Application completed in {end - start:.2f} seconds.")

    return end - start


def run_regression_test(application, new_gold_standard, test_dir):
    # Move to the application directory
    application_path = os.path.join(test_dir, application)

    # Run the simulation and move the results to the test directory
    test_time = compile_and_run_simulation(application_path)

    # Compare the result against the gold standard, if it exists
    tolerance = 1e-3
    if new_gold_standard:
        shutil.move("output.txt", "gold_output.txt")
        test_passed = True
    else:

        def extract_l2_norms(filename):
            """Extracts the l2-norm values from a given file."""
            norms = {}
            with open(filename, "r") as file:
                for line in file:
                    match = re.search(
                        r"Solution index (\d+) type \w+ l2-norm: ([\d.eE+-]+)", line
                    )
                    if match:
                        index, norm_value = int(match.group(1)), float(match.group(2))
                        norms[index] = norm_value
            return norms

        # Read values from gold standard and test output
        gold_norms = extract_l2_norms("gold_output.txt")
        test_norms = extract_l2_norms("output.txt")

        test_passed = True
        for index in gold_norms:
            if index in test_norms:
                rel_diff = abs(
                    (gold_norms[index] - test_norms[index]) / gold_norms[index]
                )
                if rel_diff >= tolerance:
                    test_passed = False
                    print(
                        f"Test failed for Solution index {index}: "
                        f"gold={gold_norms[index]}, test={test_norms[index]}, rel_diff={rel_diff}"
                    )

    # Determine test result
    test_result = (
        "New Gold Standard"
        if test_passed and new_gold_standard
        else "Pass"
        if test_passed
        else "Fail"
    )

    # Print the results to the screen
    print(
        f"Regression Test: {application}\n"
        f"Result: {test_result}\n"
        f"Time taken: {test_time}\n"
    )
    sys.stdout.flush()

    # Write the results to a file
    test_result_file = os.path.join(test_dir, "test_results.txt")
    write_to_file(
        test_result_file,
        f"Application: {application}\nResult: {test_result}\nTime: {test_time}\n \n",
    )

    return test_passed, test_time


def run_regression_tests_in_parallel(
    application_list, gold_standard_list, test_dir, n_processes=None
):
    regression_test_counter = 0
    regression_tests_passed = 0

    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = [
            executor.submit(run_regression_test, application, gold_standard, test_dir)
            for application, gold_standard in zip(application_list, gold_standard_list)
        ]

        for future in as_completed(futures):
            test_passed, test_time = future.result()
            regression_tests_passed += int(test_passed)
            regression_test_counter += 1

    return regression_tests_passed, regression_test_counter


# Initialize arg parser
parser = argparse.ArgumentParser(
    description="The maximum processes to use to run the regression tests"
)
parser.add_argument("-n", "--ntasks", type=int, default=1, help="Number of processes")
parser.add_argument(
    "-j",
    "--nthreads",
    type=int,
    default=1,
    help="Number of threads for compiling & running tests",
)
args = parser.parse_args()
n_processes = args.ntasks

# Grab current directory and the path to the test results file
pwd = os.path.dirname(os.path.realpath(__file__))
test_result_file = os.path.join(pwd, "test_results.txt")

# Grab current date
now = datetime.datetime.now()

# List of applications
applicationList = [
    "allen_cahn_explicit",
    "allen_cahn_implicit",
    "cahn_hilliard_explicit",
    "heat_equation_steady_state",
    "mechanics",
    "precipitate_explicit",
]
getNewGoldStandardList = [
    False,
    False,
    False,
    False,
    False,
    False,
]

# Grab cpu information
architecture, cpu_model, cpu_cores, cpu_max_freq, cpu_min_freq, hypervisor = (
    grab_cpu_information()
)

# Write to test results file
write_to_file(
    test_result_file,
    "--------------------------------------------------------- \n"
    "Regression test on " + now.strftime("%Y-%m-%d %H:%M") + "\n"
    f"Architecture: {architecture}\n"
    f"Model name: {cpu_model}\n"
    f"CPU(s): {cpu_cores}\n"
    f"CPU max/min MHz: {cpu_max_freq}, {cpu_min_freq}\n"
    f"Hypervisor vendor: {hypervisor}\n"
    f"Number of processes: {n_processes}\n"
    "--------------------------------------------------------- \n",
)

# Run tests in parallel
start_parallel = time.time()
regression_tests_passed, regression_test_counter = run_regression_tests_in_parallel(
    applicationList, getNewGoldStandardList, pwd, n_processes
)
end_parallel = time.time()
print(f"Total time spend on regressions tests: {end_parallel - start_parallel}")


os.chdir(pwd)
# Output the overall test results
write_to_file(
    test_result_file,
    f"Tests Passed: {regression_tests_passed}/{regression_test_counter}\n"
    "--------------------------------------------------------- \n",
)

# Print overall results
print(f"Regression Tests Passed: {regression_tests_passed}/{regression_test_counter}\n")

# Set exit code
if regression_tests_passed < regression_test_counter:
    sys.exit(1)
else:
    sys.exit(0)
