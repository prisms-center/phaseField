// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/bounding_box.h>
#include <deal.II/fe/fe_values.h>

#include <prismspf/core/cell_marker_base.h>
#include <prismspf/core/constraint_manager.h>
#include <prismspf/core/dof_manager.h>
#include <prismspf/core/system_wide.h>
#include <prismspf/core/triangulation_manager.h>

#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include "prismspf/core/field_attributes.h"

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class RefinementManager
{
public:
  /**
   * @brief Constructor. Init the flags for refinement.
   */
  explicit RefinementManager(SolverContext<dim, degree, number> &solver_context)
    : solver_context(solver_context)
    , fe_values_flags()
    , num_quad_points(SystemWide<dim, degree>::quadrature.size())
    , max_refinement(solver_context.get_user_inputs()
                       .get_spatial_discretization()
                       .get_max_refinement())
    , min_refinement(solver_context.get_user_inputs()
                       .get_spatial_discretization()
                       .get_min_refinement())
    , marker_functions()
  {
    fe_values_flags.fill(dealii::UpdateFlags::update_default);
    for (const auto &criterion : solver_context.get_user_inputs()
                                   .get_spatial_discretization()
                                   .get_refinement_criteria())
      {
        // Grab the index and field type
        const Types::Index          index = criterion.get_index();
        const FieldInfo::TensorRank rank  = solver_context.get_user_inputs()
                                             .get_variable_attributes()
                                             .at(index)
                                             .field_info.tensor_rank;

        if (criterion.get_criterion() & GridRefinement::RefinementFlags::Value)
          {
            fe_values_flags[int(rank)] |= dealii::UpdateFlags::update_values;
          }
        else if (criterion.get_criterion() & GridRefinement::RefinementFlags::Gradient)
          {
            fe_values_flags[int(rank)] |= dealii::UpdateFlags::update_gradients;
          }
      }
    // Create the FEValues
    for (const auto field_type :
         {FieldInfo::TensorRank::Scalar, FieldInfo::TensorRank::Vector})
      {
        fe_values[field_type] =
          dealii::FEValues<dim>(SystemWide<dim, degree>::fe_systems[field_type],
                                SystemWide<dim, degree>::quadrature,
                                fe_values_flags.at(int(field_type)));
      }
  }

  /**
   * @brief Destructor.
   */
  ~RefinementManager() = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so grid refiner instances aren't copied.
   */
  RefinementManager(const RefinementManager &grid_refiner) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so grid refiner instances aren't copied.
   */
  RefinementManager &
  operator=(const RefinementManager &grid_refiner) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so grid refiner instances aren't moved.
   */
  RefinementManager(RefinementManager &&grid_refiner) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so grid refiner instances aren't moved.
   */
  RefinementManager &
  operator=(RefinementManager &&grid_refiner) noexcept = delete;

  /**
   * @brief Do the adaptive refinement
   *
   * This function consists of a few steps.
   * 1. Flag the cells on the mesh according to the refinement criterion.
   * 2. Refine the grid and initialize the solution transfer object
   * 3. Redistribute the DoFs
   * 4. Transfer the solution from old to new
   * 5. Recompute and reapply the constraints (this is done in the solvers)
   * 6. Recompute invm & element volume (if applicable)
   */
  void
  do_adaptive_refinement()
  {
    // Return early if adaptive meshing is disabled
    if (!solver_context.get_user_inputs()
           .get_spatial_discretization()
           .get_has_adaptivity())
      {
        return;
      }

    TriangulationManager<dim> &triangulation_manager =
      solver_context.get_triangulation_manager();
    DofManager<dim> &dof_manager = solver_context.get_dof_manager();
    ConstraintManager<dim, degree, number> &constraint_manager =
      solver_context.get_constraint_manager();
    std::vector<FieldAttributes> &field_attributes =
      solver_context.get_field_attributes();

    // Step 1
    mark_cells_for_refinement_and_coarsening();
    bool first_iteration = true;
    while (
      dealii::Utilities::MPI::logical_or(mark_cells_for_refinement(), MPI_COMM_WORLD) ||
      first_iteration)
      {
        first_iteration = false;

        // Update ghosts of all fields.
        for (int field_index = 0; field_index < field_attributes.size(); field_index++)
          {
            solver_context.get_solution_indexer()
              .get_solution(field_index)
              .update_ghost_values();
          }

        // Step 2
        refine_grid();

        // Step 3
        triangulation_manager.reinit();
        dof_manager.reinit(triangulation_manager, field_attributes);
        constraint_manager().make_constraints(SystemWide<dim, degree>::mapping,
                                              dof_manager.get_dof_handlers());

        // Todo: reinit matrix free operators?

        // Todo: reinit solutions?

        // Step 4
        // Todo execute solution transfer

        // Step 6
        // Todo: Recompute invm & element volume
      }
  }

  void
  add_refinement_marker(std::shared_ptr<const CellMarkerBase<dim>> marker)
  {
    marker_functions.push_back(marker);
  }

  void
  clear_refinement_markers()
  {
    marker_functions.clear();
  }

  const std::vector<std::shared_ptr<const CellMarkerBase<dim>>> &
  get_refinement_markers() const
  {
    return marker_functions;
  }

private:
  /**
   * @brief Mark cells for refinement and coarsening
   */
  void
  mark_cells_for_refinement_and_coarsening()
  {
    // Create the an object for the refinement criterion at each of the quad points. This
    // will either contain the value for scalar fields, the magnitude for vector fields,
    // or the magnitude of the gradient for both of the fields.
    std::vector<number> values(num_quad_points, 0.0);

    // Clear user flags
    solver_context.get_triangulation_handler().clear_user_flags();

    // Loop over the cells provided by the triangulation
    for (const auto &cell : solver_context.get_triangulation_handler()
                              .get_triangulation()
                              .active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            // Whether we should refine the cell
            bool should_refine = false;

            // TODO (landinjm): We can probably avoid checking some of the neighboring
            // cells when coarsening them
            for (const auto &criterion : solver_context.get_user_inputs()
                                           .get_spatial_discretization()
                                           .get_refinement_criteria())
              {
                // Grab the index
                const Types::Index index = criterion.get_index();

                // Grab the field type
                const FieldInfo::TensorRank local_field_type =
                  solver_context.get_user_inputs()
                    .get_variable_attributes()
                    .at(index)
                    .field_info.tensor_rank;

                // Grab the DoFHandler iterator
                const auto dof_iterator = cell->as_dof_handler_iterator(
                  solver_context.get_dof_handler().get_dof_handler(index));

                // Reinit the cell
                fe_values.at(local_field_type).reinit(dof_iterator);

                if (criterion.get_criterion() & GridRefinement::RefinementFlags::Value)
                  {
                    if (local_field_type == FieldInfo::TensorRank::Scalar)
                      {
                        // Get the values for a scalar field
                        fe_values.at(local_field_type)
                          .get_function_values(
                            *solver_context.get_solution_handler()
                               .get_solution_vector(index, DependencyType::Normal),
                            values);
                      }
                    else
                      {
                        // Get the magnitude of the value for vector fields
                        // TODO (landinjm): Should be zeroing this out?
                        std::vector<dealii::Vector<number>> vector_values(
                          num_quad_points,
                          dealii::Vector<number>(dim));
                        fe_values.at(local_field_type)
                          .get_function_values(
                            *solver_context.get_solution_handler()
                               .get_solution_vector(index, DependencyType::Normal),
                            vector_values);
                        for (unsigned int q_point = 0; q_point < num_quad_points;
                             ++q_point)
                          {
                            values[q_point] = vector_values[q_point].l2_norm();
                          }
                      }

                    // Check if any of the quadrature points meet the refinement criterion
                    for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                      {
                        if (criterion.value_in_open_range(values[q_point]))
                          {
                            should_refine = true;
                            break;
                          }
                      }

                    // Exit if we've already determined that the cell has to be refined
                    if (should_refine)
                      {
                        break;
                      }
                  }
                if (criterion.get_criterion() & GridRefinement::RefinementFlags::Gradient)
                  {
                    if (local_field_type == FieldInfo::TensorRank::Scalar)
                      {
                        // Get the magnitude of the gradient for a scalar field
                        // TODO (landinjm): Should be zeroing this out?
                        std::vector<dealii::Tensor<1, dim, number>> scalar_gradients(
                          num_quad_points);
                        fe_values.at(local_field_type)
                          .get_function_gradients(
                            *solver_context.get_solution_handler()
                               .get_solution_vector(index, DependencyType::Normal),
                            scalar_gradients);
                        for (unsigned int q_point = 0; q_point < num_quad_points;
                             ++q_point)
                          {
                            values[q_point] = scalar_gradients[q_point].norm();
                          }
                      }
                    else
                      {
                        // TODO (landinjm): Should be zeroing this out?
                        std::vector<std::vector<dealii::Tensor<1, dim, number>>>
                          vector_gradients(num_quad_points,
                                           std::vector<dealii::Tensor<1, dim, number>>(
                                             dim));
                        fe_values.at(local_field_type)
                          .get_function_gradients(
                            *solver_context.get_solution_handler()
                               .get_solution_vector(index, DependencyType::Normal),
                            vector_gradients);
                        for (unsigned int q_point = 0; q_point < num_quad_points;
                             ++q_point)
                          {
                            dealii::Vector<number> vector_gradient_component_magnitude(
                              dim);
                            for (unsigned int dimension = 0; dimension < dim; dimension++)
                              {
                                vector_gradient_component_magnitude[dimension] =
                                  vector_gradients[q_point][dimension].norm();
                              }
                            values[q_point] =
                              vector_gradient_component_magnitude.l2_norm();
                          }
                      }

                    // Check if any of the quadrature points meet the refinement criterion
                    for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                      {
                        if (criterion.gradient_magnitude_above_threshold(values[q_point]))
                          {
                            should_refine = true;
                            break;
                          }
                      }

                    // Exit if we've already determined that the cell has to be refined
                    if (should_refine)
                      {
                        break;
                      }
                  }
              }

            Assert(cell->level() > 0,
                   dealii::ExcMessage("Cell refinement level is less than one, which "
                                      "will lead to underflow."));
            const auto cell_refinement = static_cast<unsigned int>(cell->level());
            if (should_refine && cell_refinement < max_refinement)
              {
                cell->set_user_flag();
                cell->clear_coarsen_flag();
                cell->set_refine_flag();
              }
            if (should_refine)
              {
                cell->set_user_flag();
                cell->clear_coarsen_flag();
              }
            if (!should_refine && cell_refinement > min_refinement &&
                !cell->user_flag_set())
              {
                cell->set_coarsen_flag();
              }
          }
      }
  }

  /**
   * @brief Mark cells based on function. Note: cells are only marked for refinement but
   * not coarsening.
   * @param refinement_function A function that determines if a cell should be refined.
   * @return True if any cell was marked for refinement, false otherwise.
   */
  bool
  mark_cells_for_refinement()
  {
    bool any_cell_marked = false;
    for (const auto &cell : solver_context.get_triangulation_manager()
                              .get_triangulation()
                              .active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            const auto cell_refinement = static_cast<unsigned int>(cell->level());
            if (std::any_of(
                  marker_functions.begin(),
                  marker_functions.end(),
                  [&](const std::shared_ptr<const CellMarkerBase<dim>> &marker_function)
                  {
                    return marker_function->flag(
                      *cell,
                      solver_context.get_user_inputs().get_temporal_discretization());
                  }))
              {
                cell->set_user_flag();
                cell->clear_coarsen_flag();
                if (cell_refinement < max_refinement)
                  {
                    cell->set_refine_flag();
                    any_cell_marked = true;
                  }
              }
          }
      }
    return any_cell_marked;
  }

  /**
   * @brief Refine the grid
   */
  void
  refine_grid()
  {
    // Prepare for grid refinement
    solver_context.get_triangulation_manager().prepare_for_grid_refinement();

    // Prepare the solution transfer objects
    solver_context.get_solution_handler().prepare_for_solution_transfer();

    // Execute grid refinement
    solver_context.get_triangulation_manager().execute_grid_refinement();
  }

  /**
   * @brief Grid refinement context.
   */
  SolverContext<dim, degree, number> solver_context;

  /**
   * @brief Update flags for the FEValues determined by the grid refinement
   * criterion. For now, we share one flag set for scalar fields and one for vector
   * fields.
   */
  std::array<dealii::UpdateFlags, 2> fe_values_flags;

  /**
   * @brief Finite element values for scalar and vector fields.
   */
  std::array<dealii::FEValues<dim>, 2> fe_values;

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
