// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/grid_refiner_context.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
class GridRefiner
{
public:
  /**
   * @brief Constructor.
   */
  explicit GridRefiner(GridRefinementContext<dim, degree> &grid_refinement_context)
    : grid_refinement_context(grid_refinement_context)
  {
    fe_values_flags.resize(static_cast<Types::Index>(FieldType::Vector) + 1,
                           dealii::UpdateFlags::update_default);

    for (const auto &criterion : grid_refinement_context.get_user_inputs()
                                   .get_spatial_discretization()
                                   .get_refinement_criteria())
      {
        // Grab the index and field type
        Types::Index index            = criterion.get_index();
        FieldType    local_field_type = grid_refinement_context.get_user_inputs()
                                       .get_variable_attributes()
                                       .at(index)
                                       .get_field_type();

        if (criterion.get_criterion() & GridRefinement::RefinementFlags::Value)
          {
            fe_values_flags[static_cast<Types::Index>(local_field_type)] |=
              dealii::UpdateFlags::update_values;
          }
        else if (criterion.get_criterion() & GridRefinement::RefinementFlags::Gradient)
          {
            fe_values_flags[static_cast<Types::Index>(local_field_type)] |=
              dealii::UpdateFlags::update_gradients;
          }
      }
  };

  /**
   * @brief Destructor.
   */
  ~GridRefiner() = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so grid refiner instances aren't copied.
   */
  GridRefiner(const GridRefiner &grid_refiner) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so grid refiner instances aren't copied.
   */
  GridRefiner &
  operator=(const GridRefiner &grid_refiner) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so grid refiner instances aren't moved.
   */
  GridRefiner(GridRefiner &&grid_refiner) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so grid refiner instances aren't moved.
   */
  GridRefiner &
  operator=(GridRefiner &&grid_refiner) noexcept = delete;

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
    if (!grid_refinement_context.get_user_inputs()
           .get_spatial_discretization()
           .get_has_adaptivity())
      {
        ConditionalOStreams::pout_base() << "  grid refinment disabled...\n"
                                         << std::flush;
        return;
      }

    // Step 1
    mark_cells_for_refinement_and_coarsening();

    // Step 2
    refine_grid();

    // Step 3
    grid_refinement_context.get_triangulation_handler().reinit();
    grid_refinement_context.get_dof_handler().reinit(
      grid_refinement_context.get_triangulation_handler(),
      grid_refinement_context.get_finite_element_systems(),
      grid_refinement_context.get_multigrid_info());
    grid_refinement_context.get_constraint_handler().make_constraints(
      grid_refinement_context.get_mapping(),
      grid_refinement_context.get_dof_handler().get_dof_handlers());
    if (grid_refinement_context.get_multigrid_info().has_multigrid())
      {
        const unsigned int min_level =
          grid_refinement_context.get_multigrid_info().get_mg_min_level();
        const unsigned int max_level =
          grid_refinement_context.get_multigrid_info().get_mg_max_level();
        for (unsigned int level = min_level; level <= max_level; ++level)
          {
            grid_refinement_context.get_constraint_handler().make_mg_constraints(
              grid_refinement_context.get_mapping(),
              grid_refinement_context.get_dof_handler().get_mg_dof_handlers(level),
              level);
          }
      }
    grid_refinement_context.get_matrix_free_handler().reinit(
      grid_refinement_context.get_mapping(),
      grid_refinement_context.get_dof_handler().get_dof_handlers(),
      grid_refinement_context.get_constraint_handler().get_constraints(),
      dealii::QGaussLobatto<1>(degree + 1));
    if (grid_refinement_context.get_multigrid_info().has_multigrid())
      {
        const unsigned int min_level =
          grid_refinement_context.get_multigrid_info().get_mg_min_level();
        const unsigned int max_level =
          grid_refinement_context.get_multigrid_info().get_mg_max_level();
        grid_refinement_context.get_mg_matrix_free_handler().resize(min_level, max_level);
        for (unsigned int level = min_level; level <= max_level; ++level)
          {
            grid_refinement_context.get_mg_matrix_free_handler()[level].reinit(
              grid_refinement_context.get_mapping(),
              grid_refinement_context.get_dof_handler().get_mg_dof_handlers(level),
              grid_refinement_context.get_constraint_handler().get_mg_constraints(level),
              dealii::QGaussLobatto<1>(degree + 1));
          }
      }
    grid_refinement_context.get_solution_handler().reinit(
      grid_refinement_context.get_matrix_free_handler());
    if (grid_refinement_context.get_multigrid_info().has_multigrid())
      {
        grid_refinement_context.get_solution_handler().mg_reinit(
          grid_refinement_context.get_mg_matrix_free_handler());
      }

    // Step 4
    grid_refinement_context.get_solution_handler().execute_solution_transfer();

    // Step 6
    grid_refinement_context.get_invm_handler().compute_invm();
    grid_refinement_context.get_element_volumes().compute_element_volume(
      grid_refinement_context.get_finite_element_systems().begin()->second);
  };

private:
  /**
   * @brief Mark cells for refinement and coarsening
   */
  void
  mark_cells_for_refinement_and_coarsening()
  {
    // Set quadrature rule and FEValues
    const dealii::QGaussLobatto<dim> quadrature(degree + 1);
    // TODO (landinjm): Fix for scalar and vector fields
    dealii::FESystem<dim> fe(dealii::FE_Q<dim>(dealii::QGaussLobatto<1>(degree + 1)), 1);
    dealii::FEValues<dim> fe_values(fe, quadrature, dealii::UpdateFlags::update_values);

    // Get the number of quadrature points
    const unsigned int num_quad_points = quadrature.size();

    // Get the min and max global refinements
    const unsigned int max_refinement = grid_refinement_context.get_user_inputs()
                                          .get_spatial_discretization()
                                          .get_max_refinement();
    const unsigned int min_refinement = grid_refinement_context.get_user_inputs()
                                          .get_spatial_discretization()
                                          .get_min_refinement();

    // Create the objects for the value and gradient refinements
    std::vector<double> values(num_quad_points, 0.0);

    // Clear user flags
    grid_refinement_context.get_triangulation_handler().clear_user_flags();

    // Loop over the cells provided by the triangulation
    for (const auto &cell : grid_refinement_context.get_triangulation_handler()
                              .get_triangulation()
                              .active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            // TODO (landinjm): We can probably avoid checking some of the neighboring
            // cells when coarsening them
            for (const auto &criterion : grid_refinement_context.get_user_inputs()
                                           .get_spatial_discretization()
                                           .get_refinement_criteria())
              {
                // Grab the index
                Types::Index index = criterion.get_index();

                Assert(grid_refinement_context.get_user_inputs()
                           .get_variable_attributes()
                           .at(index)
                           .get_field_type() != FieldType::Vector,
                       FeatureNotImplemented("Vector AMR"));
                Assert(criterion.get_criterion() !=
                         GridRefinement::RefinementFlags::Gradient,
                       FeatureNotImplemented("Gradient AMR criteria"));

                // Grab the DoFHandler iterator
                const auto dof_iterator = cell->as_dof_handler_iterator(
                  grid_refinement_context.get_dof_handler().get_dof_handler(index));

                // Reinit the cell
                fe_values.reinit(dof_iterator);

                // Get the values
                fe_values.get_function_values(
                  *grid_refinement_context.get_solution_handler()
                     .get_solution_vector(index, DependencyType::Normal),
                  values);

                // Check if any of the quadrature points meet the refinement criterion
                bool should_refine = false;
                for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                  {
                    if (criterion.value_in_open_range(values[q_point]))
                      {
                        should_refine = true;
                        break;
                      }
                  }

                // TODO (landinjm): This is kinda bad. Also when dealing with multiple
                // criterion we can automatically skip looking at cells if they've
                // already been flaged to be refined.
                Assert(cell->level() > 0,
                       dealii::ExcMessage("Cell refinement level is less than one, which "
                                          "will lead to underflow."));
                const auto current_cell_refinement =
                  static_cast<unsigned int>(cell->level());
                if (should_refine && current_cell_refinement < max_refinement)
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
                if (!should_refine && current_cell_refinement > min_refinement &&
                    !cell->user_flag_set())
                  {
                    cell->set_coarsen_flag();
                  }
              }
          }
      }
  }

  /**
   * @brief Refine the grid
   */
  void
  refine_grid()
  {
    // Prepare for grid refinement
    grid_refinement_context.get_triangulation_handler().prepare_for_grid_refinement();

    // Prepare the solution transfer objects
    grid_refinement_context.get_solution_handler().prepare_for_solution_transfer();

    // Execute grid refinement
    grid_refinement_context.get_triangulation_handler().execute_grid_refinement();
  }

  /**
   * @brief Grid refinement context.
   */
  GridRefinementContext<dim, degree> grid_refinement_context;

  /**
   * @brief Update flags for the FEValues object as determined by the grid refinement
   * criterion.
   */
  std::vector<dealii::UpdateFlags> fe_values_flags;
};

PRISMS_PF_END_NAMESPACE