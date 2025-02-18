import re


def parse_file(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    # Extract the number of degrees of freedom
    dof_match = re.search(r"Number of degrees of freedom:\s+(\d+)", lines[0])
    dofs = int(dof_match.group(1)) if dof_match else None

    data = {}
    for line in lines:
        match = re.match(
            r"\s*(Main init|Solve increment|Output)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)",
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


# Process multiple files
appname = "allen_cahn"
copies = "8"
filenames = [
    f"{appname}/trial_0_{copies}.txt",
    f"{appname}/trial_1_{copies}.txt",
    f"{appname}/trial_2_{copies}.txt",
]
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

# Print results
print(f"Average Degrees of Freedom: {final_dofs}")
for key, values in final_results.items():
    print(
        f"{key}: Avg Time = {values['avg_time']:.6f}, Memory = {values['memory']:.6f}"
    )
