import pyvista as pv
import numpy as np
import os
import glob
import sys

def split_mesh(input_file, output_prefix, number_of_splits):
    # Read the .vtu file
    mesh = pv.read(input_file)

    # Number of cells in each split
    cells_per_split = mesh.n_cells // number_of_splits

    split_filenames = []

    for i in range(number_of_splits):
        start = i * cells_per_split
        end = start + cells_per_split
        if i == number_of_splits - 1:
            end = mesh.n_cells

        # Extract a subset of cells
        submesh = mesh.extract_cells(np.arange(start, end))

        # Write the submesh to a new .vtu file
        split_filename = f"{output_prefix}.{i}.vtu"
        submesh.save(split_filename)
        split_filenames.append(split_filename)

    # Create the .pvtu file
    create_pvtu_file(output_prefix, split_filenames, mesh)

def create_pvtu_file(base_name, split_filenames, reference_mesh):
    with open(f"{base_name}.pvtu", 'w') as file:
        file.write('<?xml version="1.0"?>\n')
        file.write('<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian" header_type="UInt32">\n')
        file.write('  <PUnstructuredGrid GhostLevel="0">\n')
        
        # Write point data arrays
        for array in reference_mesh.point_data.keys():
            num_components = reference_mesh.point_data[array].shape[1] if len(reference_mesh.point_data[array].shape) > 1 else 1
            data_type = pv.get_array_association(reference_mesh, array)
            file.write(f'    <PDataArray type="{data_type}" Name="{array}" NumberOfComponents="{num_components}"/>\n')

        # Write piece references
        for filename in split_filenames:
            file.write(f'    <Piece Source="{filename}"/>\n')

        file.write('  </PUnstructuredGrid>\n')
        file.write('</VTKFile>\n')

def process_vtu_files(no_parts):
    # Construct the search pattern for .vtu files
    search_pattern = os.path.join('solution-*.vtu')
    
    # Find all .vtu files in the specified directory
    vtu_files = glob.glob(search_pattern)

    # Process each .vtu file sequentially
    for file in vtu_files:
        print(f"Splitting {file}")
        # Add your code here to process each file
        # For example, reading the file, performing some analysis, etc.
        filenm_noext=file[:-4]
        split_mesh(file, filenm_noext, int(no_parts))

# Main
# Call the function
try:
    process_vtu_files(sys.argv[1])
except Exception as e:
    print(f"An error occurred: {e}")
    print("Error: A single integer argument must be provided specifying the number files in which the system will be split")
    



