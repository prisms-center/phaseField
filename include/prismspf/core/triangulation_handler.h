// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria.h>

#include <prismspf/config.h>

#include <string>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
class userInputParameters;

template <int dim>
class solutionHandler;

/**
 * \brief This class handlers the generation and manipulation of triangulations.
 */
template <int dim>
class triangulationHandler
{
public:
  using Triangulation =
    std::conditional_t<dim == 1,
                       dealii::Triangulation<dim>,
                       dealii::parallel::distributed::Triangulation<dim>>;

  /**
   * \brief Constructor.
   */
  explicit triangulationHandler(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Getter function for triangulation (constant reference).
   */
  [[nodiscard]] const Triangulation &
  get_triangulation() const;

  /**
   * \brief Getter function for the multigrid triangulation (constant reference).
   */
  [[nodiscard]] const std::vector<std::shared_ptr<const dealii::Triangulation<dim>>> &
  get_mg_triangulation() const;

  /**
   * \brief Getter function for a level in the globally coarsening multigrid triangulation
   * (constant reference).
   */
  [[nodiscard]] const dealii::Triangulation<dim> &
  get_mg_triangulation(unsigned int level) const;

  /**
   * \brief Return the global maximum level of the triangulation.
   */
  [[nodiscard]] unsigned int
  get_n_global_levels() const;

  /**
   * \brief Return the maximum multigrid level.
   */
  [[nodiscard]] unsigned int
  get_mg_min_level() const;

  /**
   * \brief Return the minimum multigrid level.
   */
  [[nodiscard]] unsigned int
  get_mg_max_level() const;

  /**
   * \brief Whether multigrid has been setup.
   */
  [[nodiscard]] bool
  has_setup_multigrid() const;

  /**
   * \brief Generate mesh based on the inputs provided by the user.
   */
  void
  generate_mesh();

  /**
   * \brief Export triangulation to vtk. This is done for debugging purposes when dealing
   * with unusual meshes (e.g., circular domains).
   */
  void
  export_triangulation_as_vtk(const std::string &filename) const;

private:
  /**
   * \brief Mark the domain ids on the triangulation to get the proper mapping of
   * specified boundary conditions.
   *
   * TODO (landinjm): When the user has different meshs (or custom for that matter), we
   * should let them manually set the boundary ids.
   */
  void
  mark_boundaries() const;

  /**
   * \brief Mark certain faces of the triangulation periodic.
   */
  void
  mark_periodic();

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> *user_inputs;

  /**
   * \brief Main triangulation.
   */
  std::shared_ptr<Triangulation> triangulation;

  /**
   * \brief Collection of triangulations for each multigrid level.
   *
   * TODO (landinjm): p-multigrid
   * TODO (landinjm): Should we allow for multiple instances of multigrid? Most likely
   * not.
   */
  std::vector<std::shared_ptr<const dealii::Triangulation<dim>>> coarsened_triangulations;

  /**
   * \brief Whether we have multigrid.
   */
  bool has_multigrid = false;

  /**
   * \brief Minimum multigrid level.
   */
  unsigned int min_level = 0;

  /**
   * \brief Maximum multigrid level.
   */
  unsigned int max_level = 0;
};

PRISMS_PF_END_NAMESPACE
