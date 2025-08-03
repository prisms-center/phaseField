// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, typename number>
class MatrixFreeContainer;

template <unsigned int dim>
class MGInfo;

/**
 * @brief A little class that computes the element volume for our triangulation.
 */
template <unsigned int dim, unsigned int degree, typename number>
class ElementVolume
{
public:
  /**
   * @brief Constructor.
   */
  ElementVolume() = default;

  /**
   * @brief Initialize.
   */
  void
  initialize(
    std::shared_ptr<dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>>
      _data);

  /**
   * @brief Compute element volume for the triangulation
   */
  void
  compute_element_volume();

  /**
   * @brief Get the vector of element volumes (const reference).
   */
  [[nodiscard]] const dealii::AlignedVector<dealii::VectorizedArray<number>> &
  get_element_volumes() const;

  /**
   * @brief Get the element volume of a certain cell (const reference).
   */
  [[nodiscard]] const dealii::VectorizedArray<number> &
  get_element_volume(unsigned cell) const;

private:
  /**
   * @brief Matrix-free object.
   */
  std::shared_ptr<dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>> data;

  /**
   * @brief Vector that stores element volumes
   */
  dealii::AlignedVector<dealii::VectorizedArray<number>> element_volume;
};

/**
 * @brief A container that holds the element volumes for multigrid and normal
 * triangulations.
 */
template <unsigned int dim, unsigned int degree, typename number>
class ElementVolumeContainer
{
public:
  /**
   * @brief Constructor.
   */
  explicit ElementVolumeContainer(MGInfo<dim> &mg_info);

  /**
   * @brief Initialize the element volume container.
   */
  void
  initialize(const MatrixFreeContainer<dim, number> &matrix_free_container);

  /**
   * @brief Compute element volumes for the triangulation
   */
  void
  compute_element_volume();

  /**
   * @brief Recompute element volumes for the triangulation
   */
  void
  recompute_element_volume();

private:
  /**
   * @brief Matrix-free object handler for non-multigrid data.
   */
  ElementVolume<dim, degree, number> element_volume;

  /**
   * @brief Matrix-free object handler for multigrid data.
   */
  dealii::MGLevelObject<ElementVolume<dim, degree, float>> multigrid_element_volume;

  /**
   * @brief Min multigrid level.
   */
  unsigned int min_level = 0;

  /**
   * @brief Max multigrid level.
   */
  unsigned int max_level = 0;
};

PRISMS_PF_END_NAMESPACE
