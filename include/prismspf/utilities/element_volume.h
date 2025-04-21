// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Compute the element volume.
 */
template <unsigned int dim, unsigned int degree, typename number>
class elementVolume
{
public:
  /**
   * \brief Constructor.
   */
  elementVolume() = default;

  /**
   * \brief Initialize.
   */
  void
  initialize(
    std::shared_ptr<dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>>
      _data);

  /**
   * \brief Compute element volume for the triangulation
   */
  void
  compute_element_volume(const dealii::FESystem<dim> &fe_system);

  /**
   * \brief Get the vector of element volumes (const reference).
   */
  const dealii::AlignedVector<dealii::VectorizedArray<number>> &
  get_element_volumes() const;

private:
  /**
   * \brief Matrix-free object.
   */
  std::shared_ptr<dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>> data;

  /**
   * \brief Vector that stores element volumes
   */
  dealii::AlignedVector<dealii::VectorizedArray<number>> element_volume;
};

PRISMS_PF_END_NAMESPACE
