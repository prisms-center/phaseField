#ifndef INCLUDE_READ_VTK_H_
#define INCLUDE_READ_VTK_H_

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <boost/unordered_map.hpp>

#include <array>
#include <fstream>
#include <functional>
#include <sstream>
#include <string>
#include <utility>

// Defining a custom hash function
// namespace std
// {
//   template <int dim>
//   struct hash<dealii::Point<<dim, double>>
//   {
//     std::size_t
//     operator()(const dealii::Point<dim, double> &point) const
//     {
//       std::size_t seed = 0;
//       for (unsigned int i = 0; i < dim; ++i)
//         {
//           seed ^= std::hash<double>()(point[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
//         }
//       return seed;
//     }
//   };
// } // namespace std

/**
 * \brief A class that reads an unstructured vtk file and creates mappings between the
 * coordinates and values.
 *
 * Currently, we only support unstructured grid vtks in legacy formats. These can be
 * output from PRISMS-PF and modified with any text editor. The unstructured mesh must be
 * quad elements for 2D or hex elements for 3D.
 *
 * An arbitrary number and order of scalar and vector fields can be read-in. Note that it
 * is not necessary for the fields that you read-in to have the same variables names are
 * those in `equations.cc`.
 */
template <int dim, typename datatype>
class ReadVTK
{
  using scalar_list = std::vector<datatype>;
  using vector_list = std::vector<dealii::Vector<datatype>>;
  using scalar_map  = boost::unordered_map<dealii::Point<dim, datatype>,
                                          datatype,
                                          std::hash<dealii::Point<dim, datatype>>>;
  using vector_map  = boost::unordered_map<dealii::Point<dim, datatype>,
                                          dealii::Vector<datatype>,
                                          std::hash<dealii::Point<dim, datatype>>>;

public:
  /**
   * \brief Class constructor.
   *
   * \param filename The name of the vtk file to read-in (e.g., phi.vtk)
   */
  ReadVTK(std::string filename);

  /**
   *\brief Class destructor.
   */
  ~ReadVTK();

  /**
   * \brief Get the unique mappings between points and their scalar values.
   */
  boost::unordered_map<std::string, scalar_map>
  get_scalar_mappings();

  /**
   * \brief Get the unique mappings between points and their vector values.
   */
  boost::unordered_map<std::string, vector_map>
  get_vector_mappings();

private:
  /**
   * \brief Check for any errors in reading the line or if we are at the end of file.
   * Note that this function will only check the lines in debug mode.
   */
  void
  check_line_errors(const std::string &context);

  /**
   * \brief Check that the vtk version is correct.
   */
  void
  check_vtk_version();

  /**
   * \brief Check that the dimension of vtk matches that of the class.
   */
  void
  check_dim();

  /**
   * \brief Find the mesh type and assign the correct enum.
   */
  void
  check_mesh_type();

  /**
   * \brief Read the points from the file.
   */
  void
  read_points();

  /**
   * \brief Read the fields from the file.
   */
  void
  read_fields();

  /**
   * \brief Compute the mappings between values and points
   */
  void
  compute_mappings();

  /**
   * \brief Filename.
   */
  std::string file;

  /**
   * \brief Input stream.
   */
  std::ifstream vtk_file;

  /**
   * \brief Number of points
   */
  unsigned int n_points;

  /**
   * \brief List of points
   */
  std::vector<dealii::Point<dim, datatype>> point_list;

  /**
   * \brief List of scalar values
   */
  boost::unordered_map<std::string, scalar_list> scalar_value_list_map;

  /**
   * \brief List of vector values
   */
  boost::unordered_map<std::string, vector_list> vector_value_list_map;

  /**
   * \brief List of scalar values to point mappings
   */
  boost::unordered_map<std::string, scalar_map> scalar_value_mappings;

  /**
   * \brief List of vector values to point mappings
   */
  boost::unordered_map<std::string, vector_map> vector_value_mappings;
};

template <int dim, typename datatype>
ReadVTK<dim, datatype>::ReadVTK(std::string filename)
  : file(std::move(filename))
{
  // Open file
  vtk_file = std::ifstream(file);

  // Check that the provided file matches support verions
  check_vtk_version();

  // Check mesh type
  check_mesh_type();

  // Read points
  read_points();

  // Check dim
  check_dim();

  // Read fields
  read_fields();

  // Compute mappings
  compute_mappings();
}

template <int dim, typename datatype>
ReadVTK<dim, datatype>::~ReadVTK()
{
  vtk_file.clear();
  vtk_file.close();
}

template <int dim, typename datatype>
boost::unordered_map<std::string,
                     boost::unordered_map<dealii::Point<dim, datatype>,
                                          datatype,
                                          std::hash<dealii::Point<dim, datatype>>>>
ReadVTK<dim, datatype>::get_scalar_mappings()
{
  return scalar_value_mappings;
}

template <int dim, typename datatype>
boost::unordered_map<std::string,
                     boost::unordered_map<dealii::Point<dim, datatype>,
                                          dealii::Vector<datatype>,
                                          std::hash<dealii::Point<dim, datatype>>>>
ReadVTK<dim, datatype>::get_vector_mappings()
{
  return vector_value_mappings;
}

template <int dim, typename datatype>
void
ReadVTK<dim, datatype>::check_line_errors(const std::string &context)
{
  // Throw an error if we reach the end of file
  Assert(!vtk_file.eof(),
         dealii::ExcMessage("Error in reading vtk file. Reached end-of-file "
                            "without finding " +
                            context + " section."));

  // Throw an error if we fail to read a line
  Assert(!vtk_file.fail(),
         dealii::ExcMessage("Error in reading vtk file. Unable to read line."));
}

template <int dim, typename datatype>
void
ReadVTK<dim, datatype>::check_vtk_version()
{
  // Read file and throw an error if you can't
  AssertThrow(vtk_file.is_open(), dealii::ExcMessage("Unable to open vtk file: " + file));

  // Read the first line into a small buffer to reduce memory allocation
  std::array<char, 32> buffer {};
  vtk_file.getline(buffer.data(), buffer.size());

  // Throw an error if the version is unsupported
  AssertThrow(std::string(buffer.data()) == "# vtk DataFile Version 3.0",
              dealii::ExcMessage("vtk file is not 'vtk DataFile Version 3.0'."));
}

template <int dim, typename datatype>
void
ReadVTK<dim, datatype>::check_dim()
{
  // Loop through lines until we find CELL_TYPES.
  std::string line;
  while (std::getline(vtk_file, line))
    {
      check_line_errors("'CELL_TYPES'");

      if (line.find("CELL_TYPES") == 0)
        {
          std::getline(vtk_file, line);

          // Create string stream for line
          std::istringstream string_stream(line);

          // Find the cell type with string stream
          unsigned int cell_type = 0;
          string_stream >> cell_type;

          // Check that we have quad elements for 2D problems and hex elements for
          // 3D problems.
          AssertThrow(3 * dim + 3 == cell_type,
                      dealii::ExcMessage("The provided vtk's dimensions does match "
                                         "the one in parameters.prm."));
          break;
        }
    }
}

template <int dim, typename datatype>
void
ReadVTK<dim, datatype>::check_mesh_type()
{
  // Loop until we find DATASET
  std::string line;

  while (std::getline(vtk_file, line))
    {
      check_line_errors("'DATASET'");

      if (line.find("DATASET") == 0)
        {
          // Create string stream for line
          std::istringstream string_stream(line);

          // Find the mesh type with string stream
          std::string dataset_str;
          std::string mesh_str;
          string_stream >> dataset_str >> mesh_str;

          AssertThrow(mesh_str == "UNSTRUCTURED_GRID",
                      dealii::ExcMessage("Error in reading unstructured grid vtk file. " +
                                         mesh_str + " is currently not supported."));

          break;
        }
    }
}

template <int dim, typename datatype>
void
ReadVTK<dim, datatype>::read_points()
{
  // Loop through lines until we find POINTS.
  std::string line;

  while (std::getline(vtk_file, line))
    {
      check_line_errors("'POINTS'");

      if (line.find("POINTS") == 0)
        {
          // Create string stream for line
          std::istringstream string_stream(line);

          // Find the number of points with string stream
          std::string point_str;
          string_stream >> point_str >> n_points;

          break;
        }
    }

  // Read the points and append them to the vector point_list
  point_list.reserve(n_points);
  for (unsigned int i = 0; i < n_points; i++)
    {
      std::getline(vtk_file, line);

      // Create string stream for line and find the coordinate
      std::istringstream string_stream(line);
      datatype           x;
      datatype           y;
      datatype           z;
      string_stream >> x >> y >> z;

      switch (dim)
        {
          case 1:
            point_list.emplace_back(dealii::Point<dim, datatype>(x));
            break;
          case 2:
            point_list.emplace_back(dealii::Point<dim, datatype>(x, y));
            break;
          case 3:
            point_list.emplace_back(dealii::Point<dim, datatype>(x, y, z));
            break;
          default:
            AssertThrow(false, dealii::ExcMessage("Invalid dimension."));
            break;
        }
    }
}

template <int dim, typename datatype>
void
ReadVTK<dim, datatype>::read_fields()
{
  // Loop through the file until we reach POINT_DATA
  std::string line;
  while (std::getline(vtk_file, line))
    {
      check_line_errors("'POINT_DATA'");

      if (line.find("POINT_DATA") == 0)
        {
          break;
        }
    }

  // Identify if we have a scalar or vector field and initialize vectors to store values
  while (std::getline(vtk_file, line))
    {
      check_line_errors("'SCALARS'");
      check_line_errors("'VECTORS'");

      if (line.find("SCALARS") == 0)
        {
          // Get the variable name
          std::istringstream string_stream(line);
          std::string        scalar_str;
          std::string        var_name;
          string_stream >> scalar_str >> var_name;

          // Add vector to map
          scalar_value_list_map[var_name].reserve(n_points);

          // Skip LOOKUP_TABLE line
          std::getline(vtk_file, line);

          // Fill in vector with values. To increase robustness we loop through the lines
          // with scalar values until we encouter another SCALAR/VECTOR line or the
          // end-of-file. Then, we rewind the position by 1 so the rest of the loop
          // continues to work.
          std::streampos previous_position = vtk_file.tellg();
          while (std::getline(vtk_file, line))
            {
              check_line_errors("'SCALARS'");

              // Get the values
              std::istringstream stream(line);
              datatype           scalar_value;
              while (stream >> scalar_value) // Extract each number
                {
                  scalar_value_list_map[var_name].emplace_back(scalar_value);
                }

              if (line.find("SCALARS") == 0 || line.find("VECTORS") == 0)
                {
                  vtk_file.seekg(previous_position);
                  break;
                }
              previous_position = vtk_file.tellg();
            }
        }
      else if (line.find("VECTORS") == 0)
        {
          // Get the variable name
          std::istringstream string_stream(line);
          std::string        vector_str;
          std::string        var_name;
          string_stream >> vector_str >> var_name;

          // Add vector to map
          vector_value_list_map[var_name].reserve(n_points);

          // Fill in vector with values. To increase robustness we loop through the lines
          // with scalar values until we encouter another SCALAR/VECTOR line or the
          // end-of-file. Then, we rewind the position by 1 so the rest of the loop
          // continues to work.
          std::streampos previous_position = vtk_file.tellg();
          while (std::getline(vtk_file, line))
            {
              check_line_errors("'VECTORS'");

              // Get the values
              std::istringstream stream(line);
              datatype           x;
              datatype           y;
              datatype           z;
              string_stream >> x >> y >> z;

              dealii::Vector<datatype> vector_value;

              switch (dim)
                {
                  case 1:
                    vector_value.reinit(1);
                    vector_value[0] = x;
                    vector_value_list_map[var_name].emplace_back(vector_value);
                    break;
                  case 2:
                    vector_value.reinit(2);
                    vector_value[0] = x;
                    vector_value[1] = y;
                    vector_value_list_map[var_name].emplace_back(vector_value);
                    break;
                  case 3:
                    vector_value.reinit(3);
                    vector_value[0] = x;
                    vector_value[1] = y;
                    vector_value[2] = z;
                    vector_value_list_map[var_name].emplace_back(vector_value);
                    break;
                  default:
                    AssertThrow(false, dealii::ExcMessage("Invalid dimension."));
                    break;
                }

              if (line.find("SCALARS") == 0 || line.find("VECTORS") == 0)
                {
                  vtk_file.seekg(previous_position);
                  break;
                }
              previous_position = vtk_file.tellg();
            }
        }
    }
}

template <int dim, typename datatype>
void
ReadVTK<dim, datatype>::compute_mappings()
{
  // Create scalar value mappings by looping through each known variable
  for (const auto &entry : scalar_value_list_map)
    {
      const std::string           &var_name   = entry.first;
      const std::vector<datatype> &value_list = entry.second;

      AssertDimension(point_list.size(), value_list.size());

      scalar_map mappings;
      for (std::size_t i = 0; i < n_points; ++i)
        {
          mappings[point_list[i]] = value_list[i];
        }

      scalar_value_mappings[var_name] = std::move(mappings);
    }

  // Create vector value mappings by looping through each known variable
  for (const auto &entry : vector_value_list_map)
    {
      const std::string                           &var_name   = entry.first;
      const std::vector<dealii::Vector<datatype>> &value_list = entry.second;

      AssertDimension(point_list.size(), value_list.size());

      vector_map mappings;
      for (std::size_t i = 0; i < n_points; ++i)
        {
          mappings[point_list[i]] = value_list[i];
        }

      vector_value_mappings[var_name] = std::move(mappings);
    }
}

template class ReadVTK<1, double>;
template class ReadVTK<2, double>;
template class ReadVTK<3, double>;

#endif