// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/utilities/utilities.h>

#include <filesystem>
#include <vtkDataArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>

PRISMS_PF_BEGIN_NAMESPACE

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

  /**
   * \brief Get the names of the scalars in the vtk file.
   */
  std::vector<std::string>
  get_scalars_names();

  /**
   * \brief Get the names of the vectors in the vtk file.
   */
  std::vector<std::string>
  get_vectors_names();

  /**
   * \brief Get scalar value for a given point
   */
  double
  get_scalar_value(const dealii::Point<dim> &point, const std::string &scalar_name);

  /**
   * \brief Get vector value for a given point
   */
  dealii::Vector<double>
  get_vector_value(const dealii::Point<dim> &point, const std::string &vector_name);

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
  AssertThrow(
    output->IsHomogeneous(),
    dealii::ExcMessage(
      "The vtk file must have homogeneous cells of type VTK_HEXAHEDRON or VTK_QUAD"));

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

template <unsigned int dim>
inline std::vector<std::string>
ReadUnstructuredVTK<dim>::get_scalars_names()
{
  std::vector<std::string> scalars_names(n_scalars);
  for (unsigned int i = 0; i < n_scalars; ++i)
    {
      scalars_names[i] = reader->GetScalarsNameInFile(static_cast<int>(i));
    }
  return scalars_names;
}

template <unsigned int dim>
inline std::vector<std::string>
ReadUnstructuredVTK<dim>::get_vectors_names()
{
  std::vector<std::string> vectors_names(n_vectors);
  for (unsigned int i = 0; i < n_vectors; ++i)
    {
      vectors_names[i] = reader->GetVectorsNameInFile(static_cast<int>(i));
    }
  return vectors_names;
}

template <unsigned int dim>
inline double
ReadUnstructuredVTK<dim>::get_scalar_value(const dealii::Point<dim> &point,
                                           const std::string        &scalar_name)
{
  // Check that the scalar name is in the vtk file
  auto scalars_names = get_scalars_names();
  AssertThrow(std::find(scalars_names.begin(), scalars_names.end(), scalar_name) !=
                scalars_names.end(),
              dealii::ExcMessage(
                "The provided vtk dataset does not contain a field named " +
                scalar_name));

  // Convet the dealii point to an array
  std::array<double, 3> point_c_array = dealii_point_to_c_array<dim>(point);

  // Set the active scalar and update the reader
  reader->SetScalarsName(scalar_name.c_str());
  reader->Update();

  // Grab the output and point data
  vtkUnstructuredGrid *output     = reader->GetOutput();
  vtkPointData        *point_data = output->GetPointData();

  // Find the point id in the vtk file
  vtkIdType point_id = output->FindPoint(point_c_array.data());
  AssertThrow(point_id >= 0, dealii::ExcMessage("No matching point found in VTK grid"));

  // Get the data array
  vtkDataArray *data_array = point_data->GetArray(scalar_name.c_str());
  AssertThrow(data_array != nullptr,
              dealii::ExcMessage(std::string("Data array not found: ") + scalar_name));

  // Get the value of the scalar at the point
  return data_array->GetComponent(point_id, 0);
}

template <unsigned int dim>
inline dealii::Vector<double>
ReadUnstructuredVTK<dim>::get_vector_value(const dealii::Point<dim> &point,
                                           const std::string        &vector_name)
{
  // Check that the scalar name is in the vtk file
  auto vectors_names = get_vectors_names();
  AssertThrow(std::find(vectors_names.begin(), vectors_names.end(), vector_name) !=
                vectors_names.end(),
              dealii::ExcMessage(
                "The provided vtk dataset does not contain a field named " +
                vector_name));

  // Convet the dealii point to an array
  std::array<double, 3> point_c_array = dealii_point_to_c_array<dim>(point);

  // Set the active vector and update the reader
  reader->SetVectorsName(vector_name.c_str());
  reader->Update();

  // Grab the output and point data
  vtkUnstructuredGrid *output     = reader->GetOutput();
  vtkPointData        *point_data = output->GetPointData();

  // Find the point id in the vtk file
  vtkIdType point_id = output->FindPoint(point_c_array.data());
  AssertThrow(point_id >= 0, dealii::ExcMessage("No matching point found in VTK grid"));

  // Get the data array
  vtkDataArray *data_array = point_data->GetArray(vector_name.c_str());
  AssertThrow(data_array != nullptr,
              dealii::ExcMessage(std::string("Data array not found: ") + vector_name));

  // Get the value of the vector at the point
  dealii::Vector<double> vector_value(dim);
  for (unsigned int i = 0; i < dim; i++)
    {
      vector_value[i] = data_array->GetComponent(point_id, i);
    }
  return vector_value;
}

PRISMS_PF_END_NAMESPACE