// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/types.h>

#include <prismspf/utilities/utilities.h>
#include <prismspf/user_inputs/spatial_discretization.h>

#include <filesystem>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
class ReadFieldBase
{
public:
  /**
   * @brief Constructor
   */
  explicit ReadFieldBase(const std::string &_filename);

  /**
   * @brief Destructor
   */
  virtual ~ReadFieldBase() = default;

  /**
   * @brief Copy constructor.
   */
  ReadFieldBase(const ReadFieldBase &read_base) = delete;

  /**
   * @brief Copy assignment.
   */
  ReadFieldBase &
  operator=(const ReadFieldBase &read_base) = delete;

  /**
   * @brief Move constructor.
   */
  ReadFieldBase(ReadFieldBase &&read_base) noexcept = delete;

  /**
   * @brief Move assignment.
   */
  ReadFieldBase &
  operator=(ReadFieldBase &&read_base) noexcept = delete;

  /**
   * @brief Print the input file for debugging
   */
  virtual void
  print_file() = 0;

  /**
   * @brief Get scalar value for a given point
   */
  virtual number
  get_scalar_value(const dealii::Point<dim> &point, const std::string &scalar_name) = 0;

  /**
   * @brief Get vector value for a given point
   */
  virtual dealii::Vector<number>
  get_vector_value(const dealii::Point<dim> &point, const std::string &vector_name) = 0;

protected:
  /**
   * @brief Spatial discretization object
   */
  const SpatialDiscretization<dim> spatial_discretization;

  /**
   * @brief Filename
   */
  const std::string &filename;
};

template <unsigned int dim, typename number>
ReadFieldBase<dim, number>::ReadFieldBase(const std::string &_filename)
  : filename(_filename)
{
  // Check that the filename exists
  if (!std::filesystem::exists(filename))
    {
      AssertThrow(false, dealii::ExcMessage("File " + filename + " does not exist"));
    }
}

PRISMS_PF_END_NAMESPACE
