import re
import glob
import matplotlib.pyplot as plt


def parse_file(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    # Extract the number of degrees of freedom
    dof_match = re.search(
        r"number of degrees of freedom:\s+(\d+)", lines[5], re.IGNORECASE
    )
    if dof_match:
        dofs = int(dof_match.group(1))
    else:
        print(f"Warning: Could not find degrees of freedom in file: {filename}")
        dofs = None

    data = {}
    for line in lines:
        match = re.match(
            r"\s*(Main init|Explicit solve|Output)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
            line,
        )
        if match:
            name = match.group(1)
            avg_time = float(match.group(4))
            memory = float(match.group(6))
            if name not in data:
                data[name] = {"avg_time": [], "memory": []}
            data[name]["avg_time"].append(avg_time)
            data[name]["memory"].append(memory)

    return dofs, data


def plot(x, y, x_label, y_label, filename, plot_y_x=False):
    assert len(x) != 0
    assert len(y) != 0

    plt.clf()
    plt.figure(figsize=(6, 6))
    plt.plot(x, y, "-o", linewidth=2, label="Data")
    if plot_y_x:
        plt.plot(x, x, "--", linewidth=2, label=r"$y=x$")
        plt.legend(fontsize=16)
    plt.xlabel(x_label, fontsize=16)
    plt.ylabel(y_label, fontsize=16)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)


# Process multiple files for different copies
appname = "allen_cahn"
copy_values = ["1", "2", "4", "8", "16"]  # Example copy values, adjust as needed
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
                aggregated_data[key] = {"avg_time": [], "memory": []}
            aggregated_data[key]["avg_time"].extend(values["avg_time"])
            aggregated_data[key]["memory"].extend(values["memory"])

    # Compute averages manually
    final_dofs = sum(dof_list) / len(dof_list) if dof_list else None
    final_results = {
        key: {
            "avg_time": sum(values["avg_time"]) / len(values["avg_time"]),
            "memory": sum(values["memory"]) / len(values["memory"]),
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
            f"{key}: Avg Time = {values['avg_time']:.6f}, Memory = {values['memory']:.6f}"
        )
    print()

# PLot the results
copy_list = []
dofs_list = []
time_list = []
memory_list = []
normalized_time_list = []
normalized_memory_list = []

for copies, data in all_aggregated_data.items():
    copy_list.append(int(copies))
    dofs_list.append(data["dofs"])
    time_list.append(data["results"]["Explicit solve"]["avg_time"])
    memory_list.append(data["results"]["Explicit solve"]["memory"])

# Normalize the time and memory lists
initial_time = time_list[0] if time_list else 1
initial_memory = memory_list[0] if memory_list else 1

for time in time_list:
    normalized_time_list.append(time / initial_time)

for memory in memory_list:
    normalized_memory_list.append(memory / initial_memory)

print("Normalized Time List:", normalized_time_list)
print("Normalized Memory List:", normalized_memory_list)


plot(copy_list, time_list, "Copies", "Time (s)", "time_vs_copies.png")
plot(copy_list, memory_list, "Copies", "Memory (MB)", "memory_vs_copies.png")

plot(
    copy_list,
    normalized_time_list,
    "Copies",
    "Normalized Time",
    "normal_time_vs_copies.png",
    True,
)
plot(
    copy_list,
    normalized_memory_list,
    "Copies",
    "Normalized Memory",
    "normal_memory_vs_copies.png",
    True,
)
