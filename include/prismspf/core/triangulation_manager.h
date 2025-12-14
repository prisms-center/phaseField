// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

#include <prismspf/user_inputs/spatial_discretization.h>

#include <prismspf/config.h>

#include <string>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim>
class UserInputParameters;

template <unsigned int dim, typename number>
class SolutionHandler;

template <unsigned int dim>
class MGInfo;

/**
 * @brief This class handlers the generation and manipulation of triangulations.
 */
template <unsigned int dim>
class TriangulationManager
{
public:
  using Triangulation =
    std::conditional_t<dim == 1,
                       dealii::Triangulation<dim>,
                       dealii::parallel::distributed::Triangulation<dim>>;

  /**
   * @brief Constructor.
   */
  explicit TriangulationManager(bool _has_multigrid);

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
  reinit()
  {
    // Check that the initial global refinement matches the maximum adaptive refinement
    /* Assert(user_inputs->get_spatial_discretization().get_global_refinement() ==
             user_inputs->get_spatial_discretization().get_max_refinement(),
           dealii::ExcMessage(
             "Currently, we don't allow the initial refinement to be lower than the "
             "maximum adpative refinement level when using multigrid. This is because we "
             "have to create a sequence of coarser meshes.")); */
    if (has_multigrid)
      {
        coarsened_triangulations =
          dealii::MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
            triangulation);
        std::reverse(coarsened_triangulations.begin(), coarsened_triangulations.end());
      }
    else
      {
        coarsened_triangulations.clear();
        coarsened_triangulations.push_back(&triangulation);
      }
    // TODO (landinjm): p-multigrid
  };

  /**
   * @brief Getter function for triangulation (constant reference).
   */
  [[nodiscard]] const Triangulation &
  get_triangulation() const;

  /**
   * @brief Getter function for the multigrid triangulation (constant reference).
   */
  [[nodiscard]] const std::vector<std::shared_ptr<const dealii::Triangulation<dim>>> &
  get_mg_triangulation() const;

  /**
   * @brief Getter function for a level in the globally coarsening multigrid triangulation
   * (constant reference).
   */
  [[nodiscard]] const dealii::Triangulation<dim> &
  get_triangulation(unsigned int relative_level) const;

  /**
   * @brief Generate mesh based on the inputs provided by the user.
   */
  void
  generate_mesh(const UserInputParameters<dim> &user_inputs);

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
   * @brief Mark the domain ids on the triangulation to get the proper mapping of
   * specified boundary conditions.
   *
   * TODO (landinjm): When the user has different meshs (or custom for that matter), we
   * should let them manually set the boundary ids.
   */
  void
  mark_boundaries(const SpatialDiscretization<dim> &discretization_params) const;

  /**
   * @brief Mark certain faces of the triangulation periodic.
   */
  void
  mark_periodic(const UserInputParameters<dim> &user_inputs);

  /**
   * @brief Whether we have multigrid.
   */
  bool has_multigrid = false;

  /**
   * @brief Main triangulation.
   */
  dealii::Triangulation<dim> triangulation;

  /**
   * @brief Coarsened triangulations for multigrid.
   */
  std::vector<std::shared_ptr<const dealii::Triangulation<dim>>> coarsened_triangulations;
};

PRISMS_PF_END_NAMESPACE
