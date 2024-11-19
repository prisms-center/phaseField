import os
import subprocess
import shutil
import glob
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
    assert (
        "tests/automatic_tests" in pwd
    ), "Current directory is not within 'tests/automatic_tests'"

    # Application path assuming file structure matches GitHub repo
    app_path = pwd.replace("tests/automatic_tests", f"applications/{app_name}")

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

    # Check that the filepath exists
    does_filepath_exist(parameter_file_path)

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


def compile_and_run(app_name, new_parameter_file, test_dir):
    """Function that compile and runs the application in debug mode

    Args:
        app_name (string): Application name
        new_parameter_file (string): Parameter file
        test_dir (string): Automatic test directory

    Returns:
        String: Whether the application was able to succeed or not. The first return
        is the application name and the second is success or failure
    """
    try:
        # Navigate to tests/automatic_tests
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
            ["cmake", "."], check=True, capture_output=True, text=True
        )
        make_result = subprocess.run(
            ["make"], check=True, capture_output=True, text=True
        )

        # Run the application
        print(f"Running {app_dir} with parameter file {new_parameter_file}")
        run_result = subprocess.run(
            ["mpirun", "-n", "1", "./main-debug", "-i", new_parameter_file],
            check=True,
            capture_output=True,
            text=True,
        )

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
            compile_result.stdout + make_result.stdout + run_result.stdout,
        )

    except subprocess.CalledProcessError as exc:
        error_output = exc.stderr if exc.stderr else str(exc)
        return (app_name, f"Failed: {error_output}", "")

    except Exception as exc:
        return (app_name, f"Failed: {str(exc)}", "")


def run_tests_in_parallel(application_list):
    """Run each test in parallel

    Args:
        application_list (list): A list of application names

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
    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(compile_and_run, app_name, new_parameter_file, test_dir)
            for app_name in application_list
        ]
        for future in futures:
            app_name, status, output = future.result()
            results[app_name] = (status, output)

    return results


# Application list
application_list = [
    "CHAC_anisotropy",
    "corrosion",
    "CHAC_anisotropyRegularized",
    "corrosion_microgalvanic",
    "coupledCahnHilliardAllenCahn",
    "dendriticSolidification",
    "eshelbyInclusion",
    "fickianDiffusion",
    "grainGrowth",
    "grainGrowth_dream3d",
    "mechanics",
    "MgNd_precipitate_single_Bppp",
    "nucleationModel",
    "allenCahn",
    "nucleationModel_preferential",
    "allenCahn_conserved",
    "alloySolidification",
    "precipitateEvolution",
    "precipitateEvolution_pfunction",
    "alloySolidification_uniform",
    "spinodalDecomposition",
    "anisotropyFacet",
    "cahnHilliard",
    "CHiMaD_benchmarks/CHiMaD_benchmark1a",
    "CHiMaD_benchmarks/CHiMaD_benchmark1b",
    "CHiMaD_benchmarks/CHiMaD_benchmark1c",
    "CHiMaD_benchmarks/CHiMaD_benchmark2a",
    "CHiMaD_benchmarks/CHiMaD_benchmark3a",
    "CHiMaD_benchmarks/CHiMaD_benchmark6a",
    "CHiMaD_benchmarks/CHiMaD_benchmark7a",
]

# Run tests in parallel
results = run_tests_in_parallel(application_list)

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
