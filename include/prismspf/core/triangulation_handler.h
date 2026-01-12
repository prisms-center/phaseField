// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>

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
class TriangulationHandler
{
public:
  using Triangulation =
    std::conditional_t<dim == 1,
                       dealii::Triangulation<dim>,
                       dealii::parallel::distributed::Triangulation<dim>>;

  /**
   * @brief Constructor.
   */
  TriangulationHandler(const UserInputParameters<dim> &_user_inputs,
                       const MGInfo<dim>              &mg_info);

  /**
   * @brief Reinitialize the triangulation handler.
   *
   * This is used for AMR with multigrid so the coarsened meshes can be reinitialized.
   */
  void
  reinit()
  {
    // Create the triangulations for the coarser levels if we have at least one instance
    // of multigrid for any of the fields
    if (!has_multigrid)
      {
        return;
      }

    Assert(triangulation->n_global_levels() > 1,
           dealii::ExcMessage(
             "Multigrid preconditioners require multilevel triangulations"));

    // Check that the initial global refinement matches the maximum adaptive refinement
    Assert(user_inputs->get_spatial_discretization().get_global_refinement() ==
             user_inputs->get_spatial_discretization().get_max_refinement(),
           dealii::ExcMessage(
             "Currently, we don't allow the initial refinement to be lower than the "
             "maximum adaptive refinement level when using multigrid. This is because we "
             "have to create a sequence of coarser meshes."));

    coarsened_triangulations =
      dealii::MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
        *triangulation);

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
  get_mg_triangulation(unsigned int level) const;

  /**
   * @brief Return the maximum multigrid level.
   */
  [[nodiscard]] unsigned int
  get_mg_min_level() const;

  /**
   * @brief Return the minimum multigrid level.
   */
  [[nodiscard]] unsigned int
  get_mg_max_level() const;

  /**
   * @brief Generate mesh based on the inputs provided by the user.
   */
  void
  generate_mesh();

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
    Assert(triangulation != nullptr, dealii::ExcNotInitialized());
    triangulation->prepare_coarsening_and_refinement();
  }

  /**
   * @brief Execute grid refinement on the triangulation.
   */
  void
  execute_grid_refinement()
  {
    Assert(triangulation != nullptr, dealii::ExcNotInitialized());
    triangulation->execute_coarsening_and_refinement();
  }

  /**
   * @brief CLear all user flags.
   */
  void
  clear_user_flags()
  {
    Assert(triangulation != nullptr, dealii::ExcNotInitialized());
    triangulation->clear_user_flags();
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
  mark_boundaries() const;

  /**
   * @brief Mark certain faces of the triangulation periodic.
   */
  void
  mark_periodic();

  /**
   * @brief User-inputs.
   */
  const UserInputParameters<dim> *user_inputs;

  /**
   * @brief Main triangulation.
   */
  std::shared_ptr<Triangulation> triangulation;

  /**
   * @brief Collection of triangulations for each multigrid level.
   *
   * TODO (landinjm): p-multigrid
   * TODO (landinjm): Should we allow for multiple instances of multigrid? Most likely
   * not.
   */
  std::vector<std::shared_ptr<const dealii::Triangulation<dim>>> coarsened_triangulations;

  /**
   * @brief Whether we have multigrid.
   */
  bool has_multigrid = false;

  /**
   * @brief Minimum multigrid level.
   */
  unsigned int min_level = 0;

  /**
   * @brief Maximum multigrid level.
   */
  unsigned int max_level = 0;
};

PRISMS_PF_END_NAMESPACE
