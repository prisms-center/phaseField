// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/types.h>

#include <prismspf/utilities/utilities.h>

#include <prismspf/field_input/read_field_base.h>

#include <fstream>
#include <iostream>
#include <string>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
class ReadBinary : public ReadFieldBase<dim, number>
{
public:
  /**
   * @brief Constructor
   */
  ReadBinary(const InitialConditionFile       &_ic_file,
             const SpatialDiscretization<dim> &_spatial_discretization);

  void
  check_size(const InitialConditionFile &ic_file);

  /**
   * @brief Print the binary file to text for debugging
   */
  void
  print_file() override;

  /**
   * @brief Write a binary file for testing
   */

  static void
  write_file(const std::vector<number> &data, const InitialConditionFile &ic_file);

  /**
   * @brief Get vector value for a given point
   */
  dealii::Vector<number>
  get_vals(const std::array<dealii::types::global_dof_index, dim> &multi_index,
           const unsigned int                                      n_values);

  /**
   * @brief Get vector value for a given point
   */
  dealii::Vector<number>
  interpolate(const dealii::Point<dim> &point, const unsigned int n_pt_values);

  /**
   * @brief Get scalar value for a given point
   */
  number
  get_scalar_value(const dealii::Point<dim> &point,
                   const std::string        &scalar_name) override;

  /**
   * @brief Get vector value for a given point
   */
  dealii::Vector<number>
  get_vector_value(const dealii::Point<dim> &point,
                   const std::string        &vector_name) override;

private:
  /**
   * @brief Number of points.
   */
  dealii::types::global_dof_index n_values;

  /**
   * @brief Vector of dNdx values for each dimension.
   */
  std::array<number, dim> dNdx;

  /**
   * @brief Data array to hold the read in values.
   */
  std::vector<number> data;
};

template <unsigned int dim, typename number>
inline ReadBinary<dim, number>::ReadBinary(
  const InitialConditionFile       &_ic_file,
  const SpatialDiscretization<dim> &_spatial_discretization)
  : ReadFieldBase<dim, number>(_ic_file, _spatial_discretization)
{
  check_size(this->ic_file);
  std::ifstream dataFile(this->ic_file.filename, std::ios::binary);
  AssertThrow(dataFile,
              dealii::ExcMessage("Could not open binary file: " +
                                 this->ic_file.filename));
  data.reserve(n_values);
  const unsigned int bufsize = sizeof(number);
  char               buf[bufsize];
  number             numbuf;
  for (dealii::types::global_dof_index j = 0; j < n_values; ++j)
    {
      dataFile.read(buf, bufsize);
      memcpy(&numbuf, &buf, bufsize);
      data.push_back(numbuf);
    }
  dataFile.close();
  for (unsigned int j = 0; j < dim; ++j)
    {
      dNdx[j] = static_cast<number>(this->ic_file.n_data_points[j] - 1) /
                this->spatial_discretization.get_size()[j];
    }
}

// member functions called during construction that could belong to base or go elsewhere
template <unsigned int dim, typename number>
inline void
ReadBinary<dim, number>::check_size(const InitialConditionFile &ic_file)
{
  std::uintmax_t file_size = std::filesystem::file_size(ic_file.filename);
  // Calculate expected size: number of points * size of each point
  n_values = 1;
  for (unsigned int j = 0; j < dim; j++)
    {
      n_values *= static_cast<std::uintmax_t>(ic_file.n_data_points[j]);
    }

  std::uintmax_t expected_size = n_values * sizeof(number);
  // Make sure expected size is not zero
  AssertThrow(expected_size != 0,
              dealii::ExcMessage(
                std::string("Expected input array size is zero, check that the ") +
                std::string("number of data points in each used direction is set "
                            "correctly in the input file.")));
  // Make sure the size matches for either a scalar or vector
  AssertThrow(((file_size == expected_size) || (file_size == dim * expected_size)),
              dealii::ExcMessage(
                std::string("Expected binary file size (") +
                std::to_string(expected_size) + std::string(" bytes for scalar or ") +
                std::to_string(3 * expected_size) +
                std::string(" bytes for vector) does not match actual file size (") +
                std::to_string(file_size) + std::string(" bytes).")));
}

template <unsigned int dim, typename number>
inline void
ReadBinary<dim, number>::write_file(const std::vector<number>  &data,
                                    const InitialConditionFile &ic_file)
{
  std::ofstream dataFile(ic_file.filename, std::ios::binary);
  AssertThrow(dataFile,
              dealii::ExcMessage("Could not open binary file: " + ic_file.filename));
  const unsigned int bufsize = sizeof(number);
  char               buf[bufsize];
  for (dealii::types::global_dof_index j = 0; j < data.size(); ++j)
    {
      memcpy(&buf, &data[j], bufsize);
      dataFile.write(buf, bufsize);
    }
  dataFile.close();
}

// function to get values at specified index
// encapsulates the assumed storage order of the input data (fortran-style)
template <unsigned int dim, typename number>
inline dealii::Vector<number>
ReadBinary<dim, number>::get_vals(
  const std::array<dealii::types::global_dof_index, dim> &multi_index,
  const unsigned int                                      n_components)
{
  dealii::Vector<number> value(n_components);
  // index of the scalar or first component of vector
  dealii::types::global_dof_index flat_index    = 0;
  dealii::types::global_dof_index subspace_size = 1;
  for (unsigned int j = 0; j < dim; ++j)
    {
      flat_index += multi_index[j] * subspace_size;
      subspace_size *= this->ic_file.n_data_points[j];
    }
  flat_index *= n_components;
  for (unsigned int component = 0; component < n_components; ++component)
    {
      value[component] = this->data[flat_index + component];
    }
  return value;
}

template <unsigned int dim, typename number>
inline dealii::Vector<number>
ReadBinary<dim, number>::interpolate(const dealii::Point<dim> &point,
                                     const unsigned int        n_pt_values)
{
  std::vector<std::array<uint, 2>> indices(dim,
                                           std::array<uint, 2> {
                                             {0u, 0u}
  }); // lower and upper indices in
                                                         // each dimension
  std::vector<std::array<number, 2>> weights(
    dim,
    std::array<number, 2> {
      {number(0.0), number(0.0)}
  }); // lower and upper weights in each dimension
  for (unsigned int j = 0; j < dim; ++j)
    {
      indices[j][0] = (dealii::types::global_dof_index) std::floor(dNdx[j] * point[j]);
      AssertThrow(indices[j][0] >= 0,
                  dealii::ExcMessage(
                    "negative index encountered in ReadBinary::interpolate"));
      AssertThrow(indices[j][0] < this->ic_file.n_data_points[j],
                  dealii::ExcMessage(
                    "index out of bounds encountered in ReadBinary::interpolate"));
      if (indices[j][0] == this->ic_file.n_data_points[j] - 1)
        {
          weights[j][0] = 1.0;
          weights[j][1] = 0.0;
          indices[j][1] = indices[j][0];
        }
      else
        {
          indices[j][1] = indices[j][0] + 1;
          // Calculate the fractional position between grid points
          number pt_lower = static_cast<number>(indices[j][0]) / dNdx[j];
          weights[j][1]   = (point[j] - pt_lower) * dNdx[j]; // 1/dNdx[j] = grid spacing;
          weights[j][0]   = 1.0 - weights[j][1];
        }
    }
  dealii::Vector<number> pt_values(n_pt_values);
  pt_values = 0.0;
  // 1 << dim is 2^dim, i.e. number of corners of the cell
  for (unsigned int k = 0; k < (1 << dim); ++k)
    {
      std::array<dealii::types::global_dof_index, dim> pt_indices;
      number                                           pt_wgt = 1.0;
      for (unsigned int j = 0; j < dim; ++j)
        {
          // bit j of k tells us whether to use lower (0)
          // or upper (1) index/weight in dimension j
          bool bit      = (k >> j & 1);
          pt_indices[j] = indices[j][bit];
          pt_wgt *= weights[j][bit];
        }
      dealii::Vector<number> values_temp = this->get_vals(pt_indices, n_pt_values);
      for (unsigned int n_pt_index = 0; n_pt_index < n_pt_values; ++n_pt_index)
        {
          pt_values[n_pt_index] += pt_wgt * values_temp[n_pt_index];
        }
    }
  return pt_values;
}

template <unsigned int dim, typename number>
inline number
ReadBinary<dim, number>::get_scalar_value(const dealii::Point<dim>           &point,
                                          [[maybe_unused]] const std::string &scalar_name)
{
  return interpolate(point, 1)[0];
}

template <unsigned int dim, typename number>
inline dealii::Vector<number>
ReadBinary<dim, number>::get_vector_value(const dealii::Point<dim>           &point,
                                          [[maybe_unused]] const std::string &vector_name)
{
  return interpolate(point, dim);
}

template <unsigned int dim, typename number>
inline void
ReadBinary<dim, number>::print_file()
{
  for (dealii::types::global_dof_index j = 0; j < n_values; ++j)
    {
      ConditionalOStreams::pout_summary() << this->data[j] << "\n";
    }
}

PRISMS_PF_END_NAMESPACE
