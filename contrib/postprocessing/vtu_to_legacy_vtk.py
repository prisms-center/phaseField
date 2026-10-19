#!/usr/bin/env python3
"""Convert a .vtu/.pvtu file to legacy ASCII .vtk, writing each named
point-data array as its own SCALARS/VECTORS block.

This is deliberately NOT a thin wrapper around a generic VTK writer
(e.g. `mesh.save("out.vtk")`): when a .vtu has several named point-data
arrays and none is marked as the single "active" scalars array, VTK's
generic legacy writer dumps them all into one combined "FIELD FieldData"
block instead of individual SCALARS blocks. PRISMS-PF's own VTK reader
(include/prismspf/field_input/read_vtk.h) only recognizes the SCALARS/
VECTORS form (the same form deal.II's own write_vtk() produces), so a
FIELD-only point-data block is silently invisible to it. Writing SCALARS/
VECTORS blocks explicitly here keeps the converted file readable by
PRISMS-PF's "read initial conditions from file" / "read time from file".

Global FIELD data (e.g. PRISMS-PF's TIME/CYCLE output metadata) is
preserved at the top of the file, unaffected by the above.

Requires pyvista.
"""

import argparse

import numpy as np
import pyvista as pv


def write_array(f, values, fmt="%.17g", per_line=9):
    values = np.asarray(values).reshape(-1)
    for i in range(0, len(values), per_line):
        f.write(" ".join(fmt % v for v in values[i : i + per_line]))
        f.write("\n")


def convert(input_file, output_file):
    mesh = pv.read(input_file)
    if not mesh.IsA("vtkUnstructuredGrid"):
        mesh = mesh.cast_to_unstructured_grid()

    points = np.asarray(mesh.points)
    cells = np.asarray(mesh.cells)  # flat legacy connectivity: [n, p0..pn-1, n, ...]
    celltypes = np.asarray(mesh.celltypes)
    n_points = mesh.n_points
    n_cells = mesh.n_cells

    field_data = {
        name: np.asarray(mesh.field_data[name]).reshape(-1)
        for name in mesh.field_data.keys()
    }
    point_data = {
        name: np.asarray(mesh.point_data[name]) for name in mesh.point_data.keys()
    }

    with open(output_file, "w") as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write(f"Converted from {input_file}\n")
        f.write("ASCII\n")
        f.write("DATASET UNSTRUCTURED_GRID\n")

        if field_data:
            f.write(f"FIELD FieldData {len(field_data)}\n")
            for name, values in field_data.items():
                f.write(f"{name} 1 {len(values)} double\n")
                write_array(f, values)

        f.write(f"POINTS {n_points} double\n")
        write_array(f, points, per_line=3)

        f.write(f"CELLS {n_cells} {len(cells)}\n")
        write_array(f, cells, fmt="%d", per_line=9)

        f.write(f"CELL_TYPES {n_cells}\n")
        write_array(f, celltypes, fmt="%d", per_line=9)

        if point_data:
            f.write(f"POINT_DATA {n_points}\n")
            for name, values in point_data.items():
                if values.ndim == 1 or values.shape[1] == 1:
                    f.write(f"SCALARS {name} double 1\n")
                    f.write("LOOKUP_TABLE default\n")
                    write_array(f, values)
                elif values.shape[1] == 3:
                    f.write(f"VECTORS {name} double\n")
                    write_array(f, values, per_line=3)
                else:
                    raise ValueError(
                        f"Field '{name}' has {values.shape[1]} components; only 1 "
                        "(scalar) or 3 (vector) are supported for legacy VTK point "
                        "data."
                    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Input .vtu/.pvtu file")
    parser.add_argument("output", help="Output legacy .vtk file")
    args = parser.parse_args()
    convert(args.input, args.output)
