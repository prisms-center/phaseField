// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

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
  explicit TriangulationManager(bool _has_multigrid = false);

  explicit TriangulationManager(const SpatialDiscretization<dim> &discretization_params,
                                bool                              _has_multigrid = false)
    : TriangulationManager(_has_multigrid)
  {
    generate_mesh(discretization_params);
    reinit();
  }

  /**
   * @brief Set whether multigrid triangulations will be generated.
   */
  void
  set_using_multigrid(bool _has_multigrid)
  {
    has_multigrid = _has_multigrid;
  }

  /**
   * @brief Reinitialize the triangulation handler.
   * This is used for AMR with multigrid so the coarsened meshes can be reinitialized.
   */
  void
  reinit();

  /**
   * @brief Getter function for triangulation (constant reference).
   */
  [[nodiscard]] const Triangulation<dim> &
  get_triangulation() const;

  /**
   * @brief Getter function for a level in the globally coarsening multigrid triangulation
   * (constant reference).
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
   * @brief Whether we have multigrid.
   */
  bool has_multigrid = false;

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
};

PRISMS_PF_END_NAMESPACE
