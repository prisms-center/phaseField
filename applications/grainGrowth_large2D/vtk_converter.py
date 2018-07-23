# Script to convert a rectilinear grid vtk file to an unstructured grid vtk file

shift_x = 0.125
shift_y = 0.125
shift_z = 0

rl_grid_filename = 'resampled_2048x2048_738g.vtk'
target_variable = 'FeatureIds'

us_grid_filename = "resampled_2048x2048_738g_us.vtk"

# Parse the rectilinear grid vtk file
f = open(rl_grid_filename)

in_x_coords = False
in_y_coords = False
in_z_coords = False

in_data_list = False

skip_line = False

num_x_coords = 0
num_y_coords = 0
num_z_coords = 0

x_coords = []
y_coords = []
z_coords = []

data = []

for line in f:
    if not skip_line:

        # First make sure line isn't a comment or blank line
        stripped_line = line.strip()
        if len(stripped_line) < 1 or stripped_line[0] is "#":
            continue

        split_line = stripped_line.split()

        if split_line[0] == "DIMENSIONS":
            num_x_coords = int(split_line[1])
            num_y_coords = int(split_line[2])
            num_z_coords = int(split_line[3])

        if split_line[0] == "X_COORDINATES":
            in_x_coords = True
            continue

        if split_line[0] == "Y_COORDINATES":
            in_y_coords = True
            continue

        if split_line[0] == "Z_COORDINATES":
            in_z_coords = True
            continue

        if split_line[0] == target_variable:
            in_data_list = True
            continue

        if split_line[0] == "SCALARS" and split_line[1] == target_variable:
            in_data_list = True
            skip_line = True
            continue

        if in_x_coords:
            if not in_y_coords:
                x_coords = x_coords + split_line
            else:
                in_x_coords = False

        if in_y_coords:
            if not in_z_coords:
                y_coords = y_coords + split_line
            else:
                in_y_coords = False

        if in_z_coords:
            if split_line[0] != "POINT_DATA" and split_line[0] != "CELL_DATA":
                z_coords = z_coords + split_line
            else:
                in_z_coords = False

        if in_data_list:
            if split_line[0] != "POINT_DATA" and split_line[0] != "CELL_DATA":
                data = data + split_line
            else:
                in_data_list = False
    else:
        skip_line = False

f.close()

# Shift the mesh if necessary
x_coords = [str(float(coord) + shift_x) for coord in x_coords]
y_coords = [str(float(coord) + shift_y) for coord in y_coords]
z_coords = [str(float(coord) + shift_z) for coord in z_coords]


# Write the unstructured grid vtk file
f = open(us_grid_filename, "w")
f.write('# vtk DataFile Version 2.0\n')
f.write('converted output from ' + rl_grid_filename + '\n')
f.write('ASCII\n')
f.write('DATASET UNSTRUCTURED_GRID\n')
f.write('\n')

# Write the POINTS header
num_cells = (num_x_coords - 1) * (num_y_coords - 1)
if (num_z_coords > 1):
    num_cells = num_cells * (num_z_coords - 1)

if (num_z_coords is 1):
    num_points = 4 * num_cells
else:
    num_points = 8 * num_cells

print('num cells: ', num_cells, 'num points: ', num_points, num_x_coords, num_y_coords, num_z_coords)

f.write('POINTS ' + str(num_points) + ' double\n')

# Write the list of points for each cell
x_counter = 0
y_counter = 0
z_counter = 0

if (num_z_coords == 1):
    for y_counter in range(0, num_y_coords - 1):
        for x_counter in range(0, num_x_coords - 1):
                f.write(x_coords[x_counter] + " " + y_coords[y_counter] + " " + z_coords[z_counter] + '\n')
                f.write(x_coords[x_counter + 1] + ' ' + y_coords[y_counter] + ' ' + z_coords[z_counter] + '\n')
                f.write(x_coords[x_counter] + ' ' + y_coords[y_counter + 1] + ' ' + z_coords[z_counter] + '\n')
                f.write(x_coords[x_counter + 1] + ' ' + y_coords[y_counter + 1] + ' ' + z_coords[z_counter] + '\n')
else:
    for z_counter in range(0, num_z_coords - 1):
        for y_counter in range(0, num_y_coords - 1):
            for x_counter in range(0, num_x_coords - 1):

                f.write(x_coords[x_counter] + " " + y_coords[y_counter] + " " + z_coords[z_counter] + '\n')
                f.write(x_coords[x_counter + 1] + ' ' + y_coords[y_counter] + ' ' + z_coords[z_counter] + '\n')
                f.write(x_coords[x_counter] + ' ' + y_coords[y_counter + 1] + ' ' + z_coords[z_counter] + '\n')
                f.write(x_coords[x_counter + 1] + ' ' + y_coords[y_counter + 1] + ' ' + z_coords[z_counter] + '\n')

                f.write(x_coords[x_counter] + ' ' + y_coords[y_counter] + ' ' + z_coords[z_counter + 1] + '\n')
                f.write(x_coords[x_counter + 1] + ' ' + y_coords[y_counter] + ' ' + z_coords[z_counter + 1] + '\n')
                f.write(x_coords[x_counter] + ' ' + y_coords[y_counter + 1] + ' ' + z_coords[z_counter + 1] + '\n')
                f.write(x_coords[x_counter + 1] + ' ' + y_coords[y_counter + 1] + ' ' + z_coords[z_counter + 1] + '\n')

# Write the cell info
f.write('\nCELLS ' + str(num_cells) + ' ' + str(num_points+num_cells) + '\n')

if (num_z_coords == 1):
    points_per_cell = 4
else:
    points_per_cell = 8

point_index = 0

for y_counter in range(0, num_y_coords - 1):
    for x_counter in range(0, num_x_coords - 1):

        if (num_z_coords == 1):
            f.write(str(points_per_cell) + " " + str(point_index) + " " + str(point_index + 1) + " " + str(point_index + 3) + " " + str(point_index + 2) + "\n")
            point_index = point_index + points_per_cell
        else:
            for z_counter in range(0, num_z_coords - 1):
                f.write(str(points_per_cell) + " " + str(point_index) + " " + str(point_index + 1) + " " + str(point_index + 3) + " " + str(point_index + 2) + " " + str(point_index + 4) + " " + str(point_index + 5) + " " + str(point_index + 7) + " " + str(point_index + 6) + "\n")

                point_index = point_index + points_per_cell

# Write the cell info
f.write('\nCELL_TYPES ' + str(num_cells) + '\n')
if (num_z_coords == 1):
    cell_type = 9
else:
    cell_type = 12

for cell in range(0, num_cells):
    f.write(str(cell_type) + ' ')

# Write the point data
f.write('\n\nPOINT_DATA ' + str(num_points) + '\n')
f.write('SCALARS ' + target_variable + ' double 1\n')
f.write('LOOKUP_TABLE default\n')

if (num_z_coords == 1):

    for y_counter in range(0, num_y_coords - 1):
        for x_counter in range(0, num_x_coords - 1):

            base_index = x_counter + y_counter * (num_x_coords)

            f.write(data[base_index] + " " + data[base_index + 1] + " " + data[base_index + (num_x_coords)] + " " + data[base_index + (num_x_coords) + 1] + " " )

else:
    for z_counter in range(0, num_z_coords - 1):
        for y_counter in range(0, num_y_coords - 1):
            for x_counter in range(0, num_x_coords - 1):

                base_index = x_counter + y_counter * (num_x_coords) + z_counter * (num_x_coords*num_y_coords)

                f.write(data[base_index] + " " + data[base_index + 1] + " " + data[base_index + (num_x_coords)] + " " + data[base_index + (num_x_coords) + 1] + " " + data[base_index + num_x_coords*num_y_coords] + " " + data[base_index + num_x_coords*num_y_coords + 1] + " " + data[base_index + (num_x_coords+ num_x_coords*num_y_coords)] + " " + data[base_index + (num_x_coords) + num_x_coords*num_y_coords + 1] + " ")


f.write("\n")

f.close()
