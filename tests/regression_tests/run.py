import argparse
import subprocess
import datetime
import os
import sys
import time
from dataclasses import dataclass, field
from typing import Dict
import numpy as np
import re

# Executable name. For now we will use release mode
executable_name = "main-release"

# Relative tolerance
tolerance = 1e-4

# Mapping from the application name to the absolute directory path
app_map = {}

# Global counters for ran regression tests and passed regressions tests
passed_counter = 0
ran_counter = 0


class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


def pprint(message, type="None"):
    assert isinstance(message, str), "The message must be a string"
    assert isinstance(type, str), "The message type must be a string"

    new_message = message
    if type == "Warning":
        new_message = bcolors.WARNING + message + bcolors.ENDC
    elif type == "Info":
        new_message = bcolors.OKBLUE + message + bcolors.ENDC
    elif type == "Good":
        new_message = bcolors.OKGREEN + message + bcolors.ENDC
    elif type == "Error":
        new_message = bcolors.FAIL + message + bcolors.ENDC

    print(new_message)


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


def regression_directory_valid(regression_dir: str, applications: list = None):
    if len(regression_dir) == 0:
        pprint(f"Regression directory `{regression_dir}` is empty", "Error")
        sys.exit(1)

    abs_path = regression_dir
    if not os.path.isabs(regression_dir):
        abs_path = os.path.abspath(regression_dir)

    if not os.path.exists(abs_path):
        pprint(f"Regression directory `{regression_dir}` does not exist", "Error")
        sys.exit(1)

    if applications is not None:
        for app in applications:
            app_path = os.path.join(abs_path, app)
            if not os.path.exists(app_path):
                pprint(
                    f"Regression application directory `{app}` does not exist", "Error"
                )
                sys.exit(1)

            executable_path = os.path.join(app_path, executable_name)
            if not os.path.exists(executable_path):
                pprint(
                    f"Regression application directory `{app}` does not have a release executable `{executable_name}`",
                    "Error",
                )
                sys.exit(1)

            # Add the folder to the global mapping from app name to directory
            app_map[app] = app_path


@dataclass
class SimulationData:
    # Global metadata
    n_dofs: int = 0
    n_variables: int = 0

    # Time-dependent fields
    iterations: np.ndarray = field(default_factory=lambda: np.array([]))
    field_l2_norm: Dict[int, np.ndarray] = field(default_factory=dict)
    field_integrated_value: Dict[int, np.ndarray] = field(default_factory=dict)

    def add_iteration(self, iteration: int):
        self.iterations = np.append(self.iterations, iteration)

    def add_values(
        self,
        index: int,
        l2_norm: float,
        integrated_value: float,
    ):
        if index not in self.field_l2_norm:
            self.field_l2_norm[index] = np.array([], dtype=float)
        if index not in self.field_integrated_value:
            self.field_integrated_value[index] = np.array([], dtype=float)

        self.field_l2_norm[index] = np.append(self.field_l2_norm[index], l2_norm)
        self.field_integrated_value[index] = np.append(
            self.field_integrated_value[index], integrated_value
        )

    def print(self):
        print(f"dofs: {self.n_dofs}")
        print(f"variables: {self.n_variables}")
        print(f"iterations: {self.iterations}")
        for key, value in self.field_l2_norm.items():
            print(f"Index {key} l2-norm: {value}")
        for key, value in self.field_integrated_value.items():
            print(f"Index {key} integrated value: {value}")

    def __post_init__(self):
        if not self.field_l2_norm:
            self.field_l2_norm = {0: np.array([])}
        if not self.field_integrated_value:
            self.field_integrated_value = {0: np.array([])}

    def __eq__(self, value):
        if not isinstance(value, SimulationData):
            return False

        if self.n_dofs != value.n_dofs:
            pprint(f"DoFs not equal ({self.n_dofs} != {value.n_dofs})", "Error")
            return False
        if self.n_variables != value.n_variables:
            pprint(
                f"Number of variables not equal ({self.n_variables} != {value.n_variables})",
                "Error",
            )
            return False

        if self.iterations.size != value.iterations.size or np.all(
            self.iterations != value.iterations
        ):
            pprint(
                f"Iterations are not equal ({self.iterations} != {value.iterations})",
                "Error",
            )
            return False

        for index in range(self.n_variables):
            if self.field_l2_norm[index].size != value.field_l2_norm[
                index
            ].size or not np.allclose(
                self.field_l2_norm[index], value.field_l2_norm[index], rtol=tolerance
            ):
                pprint(
                    f"For index {index} l2-norms are not close using relative tolerance {tolerance}\n{self.field_l2_norm[index]} != {value.field_l2_norm[index]}",
                    "Error",
                )
                return False

            if self.field_integrated_value[index].size != value.field_integrated_value[
                index
            ].size or not np.allclose(
                self.field_integrated_value[index],
                value.field_integrated_value[index],
                rtol=tolerance,
            ):
                pprint(
                    f"For index {index} integrated values are not close using relative tolerance {tolerance}\n{self.field_integrated_value[index]} != {value.field_integrated_value[index]}",
                    "Error",
                )
                return False

        return True


def parse_output(file: str):
    data = SimulationData()

    with open(file, "r") as f:
        for line in f:
            line = line.strip()

            if "number of degrees of freedom:" in line.lower():
                match = re.search(
                    r"number of degrees of freedom:\s*(\d+)", line, re.IGNORECASE
                )
                if match:
                    data.n_dofs = int(match.group(1))

            elif "number of variables:" in line.lower():
                match = re.search(r"number of variables:\s*(\d+)", line, re.IGNORECASE)
                if match:
                    data.n_variables = int(match.group(1))

            elif "iteration:" in line.lower():
                match = re.search(r"iteration:\s*(\d+)", line, re.IGNORECASE)
                if match:
                    iteration = int(match.group(1))
                    data.add_iteration(iteration)

                    # Now that we've found an iteration and know the number of
                    # variables, we can add a few lines to find the l2-norms
                    # and integrated values
                    for index in range(data.n_variables):
                        try:
                            next_line = next(f).strip()
                            pattern = re.compile(
                                r"Solution\s+index\s+(\d+)\s+l2-norm:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s+integrated\s+value:\s*((?:[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?\s*)+)",
                                re.IGNORECASE,
                            )
                            match = pattern.search(next_line)
                            if match:
                                index = int(match.group(1))
                                l2_norm = float(match.group(2))
                                integrated_values = [
                                    float(x) for x in match.group(3).split()
                                ]

                                data.add_values(index, l2_norm, integrated_values)
                        except StopIteration:
                            # End of file reached
                            break

    return data


def run_application(application: str, n_threads: int, new_reference: bool = False):
    pprint(f"Running {application}...", "Info")

    global ran_counter
    ran_counter += 1

    pwd = os.getcwd()
    app_dir = app_map[application]

    os.chdir(app_dir)

    start_time = time.time()

    run_result = subprocess.run(
        ["mpirun", "-n", f"{n_threads}", f"{executable_name}"],
        capture_output=True,
        text=True,
    )

    end_time = time.time()

    passed = False
    if run_result.returncode == 0:
        with open("output.txt", "w") as f:
            f.write(run_result.stdout)

        if new_reference:
            with open("gold_output.txt", "w") as f:
                f.write(run_result.stdout)

        new_data = parse_output("output.txt")
        reference_data = parse_output("gold_output.txt")

        passed = new_data == reference_data
    else:
        print(run_result.stderr)

    os.chdir(pwd)

    if passed:
        pprint(f"passed in {end_time - start_time :.2f} seconds!", "Good")

        global passed_counter
        passed_counter += 1
    else:
        pprint("failed!", "Error")


# Initialize the argument parser and add some command line parameters
parser = argparse.ArgumentParser()
parser.add_argument(
    "-j",
    "--nthreads",
    type=int,
    default=1,
    help="Number of threads for compiling & running tests",
)
parser.add_argument(
    "-d",
    "--dir",
    type=str,
    default="install/tests/regression_tests",
    help="The directory where the regression tests are located",
)
parser.add_argument("applications", nargs="+", help="List of applications to process")
args = parser.parse_args()
n_threads = args.nthreads
apps = args.applications
dir = args.dir

# Grab cpu information
architecture, cpu_model, cpu_cores, cpu_max_freq, cpu_min_freq, hypervisor = (
    grab_cpu_information()
)

# Grab the date
date = datetime.datetime.now()

# Check that the applications exist and are compiled given the regression directory
regression_directory_valid(dir, apps)

print(90 * "=")
pprint(f"Running regression tests with {n_threads} threads per process.", "Info")
print(f"Date: {date}")
print(f"Regression Directory: {dir}")
pprint(f"Applications: {apps}", "Info")
print(90 * "-")
print(f"Machine Info")
print(f"  CPU arch: {architecture}")
print(f"  CPU model: {cpu_model}")
print(f"  CPU cores: {cpu_cores}")
print(f"  CPU min/max freq (MHz): {cpu_min_freq} / {cpu_max_freq}")
print(f"  Hypervisor: {hypervisor}")
print(90 * "-")

for app in apps:
    run_application(app, n_threads)
if passed_counter != ran_counter:
    sys.exit(1)

print(90 * "=")
