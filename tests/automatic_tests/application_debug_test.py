import os
import subprocess
import shutil
import glob
import argparse
from concurrent.futures import ProcessPoolExecutor


def does_filepath_exist(filepath):
    """Function to check that a certain filepath exists

    Args:
        filepath (string): The absolute filepath

    Raises:
        FileNotFoundError: An error if the filepath does not exist

    Returns:
        Bool: Whether the filepath exists
    """
    if os.path.exists(filepath):
        return True
    else:
        raise FileNotFoundError(f"{filepath} does not exist")


def get_application_path(app_name):
    """Function that returns that absolute path for the specified application

    Args:
        app_name (string): The application name

    Returns:
        String: Absolute application path
    """
    # Current directory
    pwd = os.getcwd()

    # Check that we're in the automatic test directory
    assert "automatic_tests" in pwd, "Current directory is not within 'automatic_tests'"

    # Application path - go up two directories from automatic_tests to get to the root
    root_dir = os.path.dirname(os.path.dirname(pwd))
    app_path = os.path.join(root_dir, "applications", app_name)

    return app_path


def set_timestep(number, app_dir, new_parameter_file):
    """Function that write create a copy of the of the parameters.prm file
       in the given application directory and set the number of timesteps

    Args:
        number (int): Number of timesteps
        app_dir (string): Absolute filepath for application directory
        new_parameter_file (string): Name of new parameters.prm
    """
    # Absolute path of the parameters file
    parameter_file_path = os.path.join(app_dir, "parameters.prm")

    # Check if the file exists, if not return early
    if not os.path.exists(parameter_file_path):
        print(f"Warning: {parameter_file_path} does not exist, skipping parameter modification")
        return

    # Make a copy of the original parameters file
    new_parameter_file_path = os.path.join(app_dir, new_parameter_file)
    shutil.copy(parameter_file_path, new_parameter_file_path)
    print(f"Copied {parameter_file_path} to {new_parameter_file_path}.")

    # Set the new number of timesteps, making sure to remove simulation end time if there
    with open(new_parameter_file_path, "r") as file:
        lines = file.readlines()

    set_n_timesteps = False
    with open(new_parameter_file_path, "w") as file:
        for line in lines:
            if "set Number of time steps" in line and not set_n_timesteps:
                file.write(f"set Number of time steps = {number}\n")
                set_n_timesteps = True
            elif "set Simulation end time" in line and not set_n_timesteps:
                file.write(f"set Number of time steps = {number}\n")
                set_n_timesteps = True
            elif "set Simulation end time" in line and set_n_timesteps:
                file.write("\n")
            else:
                file.write(line)
    print(
        f"Updated {new_parameter_file} in {app_dir} to set Number of time steps to {number}."
    )


def compile_and_run(app_name, new_parameter_file, test_dir, n_threads=1, run_application=True):
    """Function that compile and runs the application in debug mode

    Args:
        app_name (string): Application name
        new_parameter_file (string): Parameter file
        test_dir (string): Automatic test directory
        n_threads (int): Number of threads to use
        run_application (bool): Whether to run the application after compilation

    Returns:
        String: Whether the application was able to succeed or not. The first return
        is the application name and the second is success or failure
    """
    try:
        # Navigate to automatic_tests
        os.chdir(test_dir)

        # Grab application path
        app_dir = get_application_path(app_name)

        # Navigate to application directory
        os.chdir(app_dir)
        print(f"Currently in {os.getcwd()}")

        # Compile the application in debug mode
        if os.path.exists("CMakeCache.txt"):
            os.remove("CMakeCache.txt")
        print(f"Compiling {app_dir}")
        compile_result = subprocess.run(
            ["cmake", ".", "-G", "Ninja"], check=True, capture_output=True, text=True
        )
        make_result = subprocess.run(
            ["ninja", "-j", str(n_threads)], check=True, capture_output=True, text=True
        )

        # Run the application if requested
        if run_application:
            print(f"Running {app_dir} with parameter file {new_parameter_file}")
            run_result = subprocess.run(
                ["mpirun", "-n", str(n_threads), "./main-debug", "-i", new_parameter_file],
                check=True,
                capture_output=True,
                text=True,
            )
            output = compile_result.stdout + make_result.stdout + run_result.stdout
        else:
            output = compile_result.stdout + make_result.stdout

        # Clean up
        print(f"Cleaning up {app_dir}")
        os.remove("CMakeCache.txt")
        os.remove("main")
        if os.path.exists("integratedFields.txt"):
            os.remove("integratedFields.txt")
        for solution_file in glob.glob("solution-*.vtu"):
            os.remove(solution_file)
        for restart_file in glob.glob("restart.*"):
            os.remove(restart_file)
        os.remove(new_parameter_file)

        # Return success
        return (
            app_name,
            "Success",
            output,
        )

    except subprocess.CalledProcessError as exc:
        error_output = exc.stderr if exc.stderr else str(exc)
        return (app_name, f"Failed: {error_output}", "")

    except Exception as exc:
        return (app_name, f"Failed: {str(exc)}", "")


def run_tests_in_parallel(application_list, n_processes, n_threads, run_application=True):
    """Run each test in parallel

    Args:
        application_list (list): A list of application names
        n_processes (int): Number of processes to use
        n_threads (int): Number of threads to use
        run_application (bool): Whether to run the application after compilation

    Returns:
        Dict: Dictionary that track the success/failure of each application
    """
    # Parameter file name
    new_parameter_file = "parameters_new.prm"

    # Grab currnet directory which should be the automatic test directory
    test_dir = os.getcwd()

    # Update parameter file for one timestep
    for app_name in application_list:
        # Grab application path
        app_dir = get_application_path(app_name)

        # Set timestep
        set_timestep(1, app_dir, new_parameter_file)

    # Dictionary for results
    results = {}

    # Run each task in parallel
    with ProcessPoolExecutor(max_workers=n_processes) as executor:
        futures = [
            executor.submit(compile_and_run, app_name, new_parameter_file, test_dir, n_threads, run_application)
            for app_name in application_list
        ]
        for future in futures:
            app_name, status, output = future.result()
            results[app_name] = (status, output)

    return results


# Initialize arg parser
parser = argparse.ArgumentParser()
parser.add_argument("-n", "--ntasks", type=int, default=1, help="Number of processes")
parser.add_argument(
    "-j",
    "--nthreads",
    type=int,
    default=1,
    help="Number of threads for compiling & running tests",
)
parser.add_argument(
    "--no-run",
    action="store_true",
    help="Only compile the applications without running them",
)
args = parser.parse_args()
n_processes = args.ntasks
n_threads = args.nthreads
run_application = not args.no_run

print(f"Running application debug tests with {n_processes} processes and {n_threads} threads.")
if not run_application:
    print("Only compiling applications, not running them.")

# Application list
application_list = [
    "allen_cahn_explicit",
    "allen_cahn_implicit",
]

# Run tests in parallel
results = run_tests_in_parallel(application_list, n_processes, n_threads, run_application)

# Print the results
print("\n\nCompilation and Execution Results:")
passed_count = 0

for app_name, (status, output) in results.items():
    print(f"{app_name}: {status}")
    if status == "Success":
        passed_count += 1
    else:
        print(f"Error details:\n{output}\n")

print(f"\nTotal applications passed: {passed_count} out of {len(results)}")
