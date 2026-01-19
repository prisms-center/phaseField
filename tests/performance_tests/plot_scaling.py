import glob
import re

import matplotlib.pyplot as plt


def parse_file(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    # Extract the number of degrees of freedom
    dofs = None
    for line in lines:
        dof_match = re.search(r"number of degrees of freedom:\s+(\d+)", line)
        if dof_match:
            dofs = int(dof_match.group(1))
            break
    if dofs is None:
        print(f"Warning: Could not find degrees of freedom in file: {filename}")

    # Grab the Caliper data
    data = {}
    found_caliper_output = False
    for line in lines:
        # Check for the beginning of the caliper output
        if re.match("Path.+", line):
            found_caliper_output = True

        if found_caliper_output:
            # Extract name, average time, and allocated memory
            match = re.match(
                r"\s*(\w[\w\s]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
                line,
            )
            if match:
                name = match.group(1).strip()
                avg_time = float(match.group(4))
                allocated_mb = float(match.group(6))
                if name not in data:
                    data[name] = {"avg_time": [], "allocated_mb": []}
                data[name]["avg_time"].append(avg_time)
                data[name]["allocated_mb"].append(allocated_mb)

    return dofs, data


def plot(x, y, x_label, y_label, filename, plot_y_x=False):
    assert len(x) != 0
    assert len(y) != 0

    plt.clf()
    plt.figure(figsize=(4, 4))
    plt.plot(x, y, "-o", linewidth=2, label="Data")
    if plot_y_x:
        plt.plot(x, x, "--", linewidth=2, label=r"$y=x$")
        plt.legend(fontsize=16)
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close("all")


# Process multiple files for different copies
appname = "allen_cahn"
copy_values = [
    "1",
    "2",
    "4",
    "8",
    "16",
    "32",
    "64",
]
all_aggregated_data = {}

for copies in copy_values:
    pattern = f"{appname}/trial_*_{copies}.txt"
    filenames = glob.glob(pattern)
    dof_list = []
    aggregated_data = {}

    for file in filenames:
        dofs, results = parse_file(file)
        if dofs:
            dof_list.append(dofs)

        for key, values in results.items():
            if key not in aggregated_data:
                aggregated_data[key] = {"avg_time": [], "allocated_mb": []}
            aggregated_data[key]["avg_time"].extend(values["avg_time"])
            aggregated_data[key]["allocated_mb"].extend(values["allocated_mb"])

    # Compute averages manually
    final_dofs = sum(dof_list) / len(dof_list) if dof_list else None
    final_results = {
        key: {
            "avg_time": sum(values["avg_time"]) / len(values["avg_time"]),
            "allocated_mb": sum(values["allocated_mb"]) / len(values["allocated_mb"]),
        }
        for key, values in aggregated_data.items()
    }

    all_aggregated_data[copies] = {"dofs": final_dofs, "results": final_results}

# Print results for each copy value
for copies, data in all_aggregated_data.items():
    print(f"Copies: {copies}")
    print(f"Average Degrees of Freedom: {data['dofs']}")
    for key, values in data["results"].items():
        print(
            f"{key}: Avg Time = {values['avg_time']:.6f}, Allocated MB = {values['allocated_mb']:.6f}"
        )
    print()

# Plot the results for each key
for key in all_aggregated_data[copy_values[0]]["results"].keys():
    copy_list = []
    dofs_list = []
    time_list = []
    memory_list = []
    normalized_time_list = []
    normalized_memory_list = []

    for copies, data in all_aggregated_data.items():
        copy_list.append(int(copies))
        dofs_list.append(data["dofs"])
        time_list.append(data["results"][key]["avg_time"])
        memory_list.append(data["results"][key]["allocated_mb"])

    # Normalize the time and memory lists
    initial_time = time_list[0] if time_list else 1
    initial_memory = memory_list[0] if memory_list else 1

    for time in time_list:
        normalized_time_list.append(time / initial_time)

    for memory in memory_list:
        normalized_memory_list.append(memory / initial_memory)

    print(f"Normalized Time List for {key}:", normalized_time_list)
    print(f"Normalized Memory List for {key}:", normalized_memory_list)

    plot(
        copy_list, time_list, "Copies", f"Time (s) - {key}", f"time_vs_copies_{key}.png"
    )
    plot(
        copy_list,
        memory_list,
        "Copies",
        f"Memory (MB) - {key}",
        f"memory_vs_copies_{key}.png",
    )

    plot(
        copy_list,
        normalized_time_list,
        "Copies",
        f"Normalized Time - {key}",
        f"normal_time_vs_copies_{key}.png",
        True,
    )
    plot(
        copy_list,
        normalized_memory_list,
        "Copies",
        f"Normalized Memory - {key}",
        f"normal_memory_vs_copies_{key}.png",
        True,
    )
