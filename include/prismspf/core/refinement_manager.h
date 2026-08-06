// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/bounding_box.h>
#include <deal.II/fe/fe_values.h>

#include <prismspf/core/cell_marker_base.h>
#include <prismspf/core/constraint_manager.h>
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/matrix_free_manager.h>
#include <prismspf/core/system_wide.h>
#include <prismspf/core/triangulation_manager.h>

#include <prismspf/solvers/solve_context.h>
#include <prismspf/solvers/solver_base.h>

#include <prismspf/user_inputs/spatial_discretization.h>

#include <prismspf/utilities/timer.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class RefinementManager
{
public:
  /**
   * @brief Constructor. Init the flags for refinement.
   */
  explicit RefinementManager(SolveContext<dim, degree, number> &_solve_context);

  ~RefinementManager() = default;

  RefinementManager(const RefinementManager &grid_refiner) = delete;
  RefinementManager &
  operator=(const RefinementManager &grid_refiner)             = delete;
  RefinementManager(RefinementManager &&grid_refiner) noexcept = delete;
  RefinementManager &
  operator=(RefinementManager &&grid_refiner) noexcept = delete;

  /**
   * @brief Do the adaptive refinement
   *
   * Perform a loop of flagging cells for refinement/coarsening and refining until no more
   * cells are flagged.
   */
  void
  do_adaptive_refinement(
    std::vector<std::shared_ptr<SolverBase<dim, degree, number>>> &solvers);

  void
  add_refinement_marker(std::shared_ptr<const CellMarkerBase<dim>> marker);

  void
  clear_refinement_markers();

  const std::list<std::shared_ptr<const CellMarkerBase<dim>>> &
  get_refinement_markers() const;

  /**
   * @brief Similar to `do_adaptive_refinement` but loops coarsening.
   *
   * Perform a loop of flagging cells for refinement/coarsening and refining until no more
   * cells are flagged.
   */
  void
  do_initial_refinement(
    std::vector<std::shared_ptr<SolverBase<dim, degree, number>>> &solvers);

private:
  /**
   * @brief Mark cells for refinement and coarsening
   */
  void
  mark_cells_for_refinement_and_coarsening();

  /**
   * @brief Mark cells based on function. Note: cells are only marked for refinement but
   * not coarsening.
   * @param refinement_function A function that determines if a cell should be refined.
   * @return True if any cell was marked for refinement, false otherwise.
   */
  bool
  mark_cells_for_refinement();

  /**
   * @brief Refine the grid once.
   */
  void
  refine_grid(std::vector<std::shared_ptr<SolverBase<dim, degree, number>>> &solvers);

  /**
   * @brief Grid refinement context.
   */
  SolveContext<dim, degree, number> *solve_context;

  /**
   * @brief Map of field name to field index
   *
   * @todo Can eliminate the map to something a little better during the constructor for
   * performance
   */
  std::map<std::string, Types::Index> field_indices;

  /**
   * @brief Update flags for the FEValues determined by the grid refinement
   * criterion. For now, we share one flag set for scalar fields and one for vector
   * fields.
   */
  std::array<dealii::UpdateFlags, 2> fe_values_flags;

  /**
   * @brief Number of quadrature points.
   */
  unsigned int num_quad_points = 0;

  /**
   * @brief Maximum global refinement level.
   */
  unsigned int max_refinement = 0;

  /**
   * @brief Minimum global refinement level.
   */
  unsigned int min_refinement = 0;

  /**
   * @brief Marker functions.
   */
  std::list<std::shared_ptr<const CellMarkerBase<dim>>> marker_functions;
};

PRISMS_PF_END_NAMESPACE
