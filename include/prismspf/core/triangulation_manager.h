// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria.h>

#include <prismspf/user_inputs/constraint_parameters.h>
#include <prismspf/user_inputs/spatial_discretization.h>

#include <prismspf/config.h>

#include <memory>
#include <string>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class handlers the generation and manipulation of triangulations.
 */
template <unsigned int dim>
class TriangulationManager
{
public:
  /**
   * @brief Constructor.
   */
  explicit TriangulationManager();

  /**
   * @brief Initialize the multigrid triangulations.
   */
  void
  init_mg();

  /**
   * @brief Remove the multigrid triangulations.
   */
  void
  clear_mg();

  /**
   * @brief Check if the multigrid triangulations are initialized.
   */
  [[nodiscard]] bool
  has_mg() const;

  /**
   * @brief Getter function for triangulation.
   */
  [[nodiscard]] const Triangulation<dim> &
  get_triangulation() const;

  /**
   * @brief Getter function for a relative level in the globally coarsening multigrid
   * triangulation.
   */
  [[nodiscard]] const dealii::Triangulation<dim> &
  get_triangulation(unsigned int relative_level) const;

  /**
   * @brief Get the vector of periodic face pairs
   */
  [[nodiscard]] const std::vector<dealii::GridTools::PeriodicFacePair<
    typename dealii::Triangulation<dim>::cell_iterator>> &
  get_periodic_face_pairs() const;

  /**
   * @brief Get the volume of the triangulation.
   */
  [[nodiscard]] double
  get_volume() const;

  /**
   * @brief Generate mesh based on the inputs provided by the user.
   */
  void
  generate_mesh(const SpatialDiscretization<dim> &discretization_params);

  /**
   * @brief Export triangulation to vtk. This is done for debugging purposes when dealing
   * with unusual meshes (e.g., circular domains).
   */
  void
  export_triangulation_as_vtk(const std::string &filename) const;

  /**
   * @brief Prepare the triangulation for grid refinement.
   */
  void
  prepare_for_grid_refinement()
  {
    triangulation.prepare_coarsening_and_refinement();
  }

  /**
   * @brief Execute grid refinement on the triangulation.
   */
  void
  execute_grid_refinement()
  {
    triangulation.execute_coarsening_and_refinement();
  }

  /**
   * @brief Clear all user flags.
   */
  void
  clear_user_flags()
  {
    triangulation.clear_user_flags();
  }

private:
  /**
   * @brief Main triangulation.
   */
  Triangulation<dim> triangulation;

  /**
   * @brief Coarsened triangulations for multigrid.
   */
  std::vector<std::shared_ptr<const dealii::Triangulation<dim>>> coarsened_triangulations;

  /**
   * @brief Periodic face pairs on the coarsest mesh if they exist
   */
  std::vector<dealii::GridTools::PeriodicFacePair<
    typename dealii::Triangulation<dim>::cell_iterator>>
    periodicity_vector;

  /**
   * @brief Volume of the triangulation.
   */
  double volume = 0;
};

PRISMS_PF_END_NAMESPACE
