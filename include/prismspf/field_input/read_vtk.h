// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <filesystem>
#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

template <unsigned int dim>
class ReadUnstructuredVTK
{
public:
  /**
   * \brief Constructor
   */
  explicit ReadUnstructuredVTK(const std::string &filename);

  /**
   * \brief Destructor
   */
  ~ReadUnstructuredVTK() = default;

  /**
   * \brief Get the vtk output
   */
  vtkUnstructuredGrid *
  get_output();

  /**
   * \brief Get the number of points
   */
  [[nodiscard]] dealii::types::global_dof_index
  get_n_points() const;

  /**
   * \brief Get the number of cells
   */
  [[nodiscard]] dealii::types::global_dof_index
  get_n_cells() const;

  /**
   * \brief Print the vtk file for debugging
   */
  void
  print_vtk_file();

private:
  /**
   * \brief Reader for the vtk file
   */
  vtkNew<vtkUnstructuredGridReader> reader;

  /**
   * \brief Number of points.
   */
  dealii::types::global_dof_index n_points;

  /**
   * \brief Number of cells.
   */
  dealii::types::global_dof_index n_cells;

  /**
   * \brief Number of scalars in file.
   */
  unsigned int n_scalars;

  /**
   * \brief Number of vectors in file.
   */
  unsigned int n_vectors;
};

template <unsigned int dim>
ReadUnstructuredVTK<dim>::ReadUnstructuredVTK(const std::string &filename)
{
  // Check that the filename exists
  if (!std::filesystem::exists(filename))
    {
      AssertThrow(false, dealii::ExcMessage("File " + filename + " does not exist"));
    }

  // Create a reader for the vtk file and update it
  // vtkNew is a smart pointer so we don't need to manage it with delete
  reader = vtkNew<vtkUnstructuredGridReader>();
  reader->SetFileName(filename.c_str());
  reader->Update();

  // Check that the file is an unstructured grid
  AssertThrow(reader->IsFileUnstructuredGrid(),
              dealii::ExcMessage("The vtk file must be an unstructured grid"));

  // Check that we only have one cell type
  auto output = reader->GetOutput();
  AssertThrow(output->IsHomogeneous(),
              dealii::ExcMessage(
                "The vtk file must have homogeneous cells of type VTK_HEXAHEDRON"));

  // Check that the cells are hexahedra or quads
  if constexpr (dim == 3)
    {
      AssertThrow(output->GetCellType(0) == VTK_HEXAHEDRON,
                  dealii::ExcMessage(
                    "For 3D meshes, the cells must be of type VTK_HEXAHEDRON "));
    }
  else if constexpr (dim == 2)
    {
      AssertThrow(output->GetCellType(0) == VTK_QUAD,
                  dealii::ExcMessage(
                    "For 2D meshes, the cells must be of type VTK_QUAD"));
    }
  else
    {
      AssertThrow(false,
                  dealii::ExcMessage("File read-in is not supported for 1D meshes"));
    }

  // Get the number of points and cells. We first fill the variables in the same type as
  // the VTK return type so we can check for types mismatches with deal.II
  vtkIdType n_points_vtk = output->GetNumberOfPoints();
  vtkIdType n_cells_vtk  = output->GetNumberOfCells();

  // Check that the number of points and cells are not too large
  AssertThrow(n_points_vtk < std::numeric_limits<dealii::types::global_dof_index>::max(),
              dealii::ExcMessage(
                "The number of points being read-in from the vtk file is too large. Try "
                "recompiling deal.II with 64-bit indices."));
  AssertThrow(n_cells_vtk < std::numeric_limits<dealii::types::global_dof_index>::max(),
              dealii::ExcMessage(
                "The number of cells being read-in from the vtk file is too large. Try "
                "recompiling deal.II with 64-bit indices."));

  // Convert the number of points and cells to the dealii type
  n_points = n_points_vtk;
  n_cells  = n_cells_vtk;

  // Get the number of scalars and vectors
  n_scalars = reader->GetNumberOfScalarsInFile();
  n_vectors = reader->GetNumberOfVectorsInFile();
}

template <unsigned int dim>
inline vtkUnstructuredGrid *
ReadUnstructuredVTK<dim>::get_output()
{
  return reader->GetOutput();
}

template <unsigned int dim>
inline dealii::types::global_dof_index
ReadUnstructuredVTK<dim>::get_n_points() const
{
  return n_points;
}

template <unsigned int dim>
inline dealii::types::global_dof_index
ReadUnstructuredVTK<dim>::get_n_cells() const
{
  return n_cells;
}

template <unsigned int dim>
inline void
ReadUnstructuredVTK<dim>::print_vtk_file()
{
  // TODO (landinjm): Should we print only on rank 0?.
  reader->GetOutput()->PrintSelf(std::cout, vtkIndent());
}