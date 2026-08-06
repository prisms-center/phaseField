#include <prismspf/core/refinement_manager.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
RefinementManager<dim, degree, number>::RefinementManager(
  SolveContext<dim, degree, number> &_solve_context)
  : solve_context(&_solve_context)
  , field_indices(field_index_map(_solve_context.get_field_attributes()))
  , fe_values_flags()
  , num_quad_points(SystemWide<dim, degree>::quadrature.size())
  , max_refinement(solve_context->get_user_inputs().spatial_discretization.max_refinement)
  , min_refinement(solve_context->get_user_inputs().spatial_discretization.min_refinement)
  , marker_functions()
{
  fe_values_flags.fill(dealii::UpdateFlags::update_default);
  std::map<std::string, Types::Index> field_indices =
    field_index_map(solve_context->get_field_attributes());
  for (const auto &[name, field_criterion] :
       solve_context->get_user_inputs().spatial_discretization.refinement_criteria)
    {
      // Grab the index and field type
      const unsigned   field_index = field_indices.at(name);
      const TensorRank rank =
        solve_context->get_field_attributes().at(field_index).field_type;

      if (field_criterion.criterion & RefinementFlags::Value)
        {
          fe_values_flags[int(rank)] |= dealii::UpdateFlags::update_values;
        }
      else if (field_criterion.criterion & RefinementFlags::Gradient)
        {
          fe_values_flags[int(rank)] |= dealii::UpdateFlags::update_gradients;
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
RefinementManager<dim, degree, number>::do_adaptive_refinement(
  std::vector<std::shared_ptr<SolverBase<dim, degree, number>>> &solvers)
{
  // Return early if adaptive meshing is disabled
  if (!solve_context->get_user_inputs().spatial_discretization.has_adaptivity)
    {
      return;
    }

  Timer::Scope timer_scope("Refine Grid");

  // Mark cells and refine until no more cells are marked
  mark_cells_for_refinement_and_coarsening();
  bool first_iteration = true;
  while (
    dealii::Utilities::MPI::logical_or(mark_cells_for_refinement(), MPI_COMM_WORLD) ||
    first_iteration)
    {
      first_iteration = false;
      refine_grid(solvers);
    }

  // Update anything affected by the grid change:
  // Recalculate InvM since the grid has changed
  solve_context->get_invm_manager().reinit(solve_context->get_matrix_free_manager());
  solve_context->get_invm_manager().compute_invm();
  // Update the ghosts
  for (auto &solver : solvers)
    {
      solver->update_ghosts();
    }

  Logger::instance()
    << LogFormatter::debug(
         "[Iteration " +
         std::to_string(solve_context->get_simulation_timer().get_increment()) +
         "] DoFs " + std::to_string(solve_context->get_dof_manager().get_total_dofs()))
    << std::endl;
}

template <unsigned int dim, unsigned int degree, typename number>
void
RefinementManager<dim, degree, number>::add_refinement_marker(
  std::shared_ptr<const CellMarkerBase<dim>> marker)
{
  marker_functions.push_back(marker);
}

template <unsigned int dim, unsigned int degree, typename number>
void
RefinementManager<dim, degree, number>::clear_refinement_markers()
{
  marker_functions.clear();
}

template <unsigned int dim, unsigned int degree, typename number>
const std::list<std::shared_ptr<const CellMarkerBase<dim>>> &
RefinementManager<dim, degree, number>::get_refinement_markers() const
{
  return marker_functions;
}

template <unsigned int dim, unsigned int degree, typename number>
void
RefinementManager<dim, degree, number>::do_initial_refinement(
  std::vector<std::shared_ptr<SolverBase<dim, degree, number>>> &solvers)
{
  Timer::Scope timer_scope("Initial Refinement");

  Logger::instance() << "Refining grid..." << std::endl;
  {
    Logger::IndentScope indent;

    const SpatialDiscretization<dim> &space_parameters =
      solve_context->get_user_inputs().spatial_discretization;
    const DoFManager<dim, degree> &dof_manager = solve_context->get_dof_manager();

    dealii::types::global_dof_index old_dofs = dof_manager.get_total_dofs();
    dealii::types::global_dof_index new_dofs = 0;

    Logger::instance() << "Initial DoFs " << old_dofs << std::endl;

    for (unsigned int remesh_index = 0; remesh_index < (space_parameters.max_refinement -
                                                        space_parameters.min_refinement);
         remesh_index++)
      {
        Logger::instance() << "Performing grid refinement..." << std::endl;
        Logger::IndentScope indent_2;

        // Perform grid refinement
        do_adaptive_refinement(solvers);

        // Update the ghosts
        for (auto &solver : solvers)
          {
            solver->update_ghosts();
          }

        // Recalculate the total DoFs
        new_dofs = dof_manager.get_total_dofs();

        // Check for convergence
        if (old_dofs == new_dofs)
          {
            break;
          }
        old_dofs = new_dofs;
      }

    Logger::instance() << "Final DoFs " << old_dofs << std::endl;
  }
  Logger::instance() << LogFormatter::success("Refinement succeeded") << std::endl;
}

template <unsigned int dim, unsigned int degree, typename number>
void
RefinementManager<dim, degree, number>::mark_cells_for_refinement_and_coarsening()
{
  Timer::Scope timer_scope("mark_cells_for_refinement_and_coarsening()");

  // Create the an object for the refinement criterion at each of the quad points. This
  // will either contain the value for scalar fields, the magnitude for vector fields,
  // or the magnitude of the gradient for both of the fields.
  std::vector<number> values(num_quad_points, 0.0);

  // Clear user flags
  solve_context->get_triangulation_manager().clear_user_flags();

  // Loop over the cells provided by the triangulation
  for (const auto &cell : solve_context->get_triangulation_manager()
                            .get_triangulation()
                            .active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Whether we should refine the cell
          bool should_refine = false;

          // TODO (landinjm): We can probably avoid checking some of the neighboring
          // cells when coarsening them
          for (const auto &[name, field_criterion] :
               solve_context->get_user_inputs()
                 .spatial_discretization.refinement_criteria)
            {
              // Grab the index
              const unsigned int index = field_indices.at(name);

              // Grab the field type
              const TensorRank local_field_type =
                solve_context->get_field_attributes().at(index).field_type;

              // Grab the DoFHandler iterator
              const auto dof_iterator = cell->as_dof_handler_iterator(
                solve_context->get_dof_manager().get_field_dof_handler(index));

              // Reinit the cell
              dealii::FEValues<dim> fe_values(
                SystemWide<dim, degree>::fe_systems[int(local_field_type)],
                SystemWide<dim, degree>::quadrature,
                fe_values_flags.at(int(local_field_type)));

              fe_values.reinit(dof_iterator);
              const auto &solution_vector =
                solve_context->get_solution_indexer().get_solution_vector(index);

              if (field_criterion.criterion & RefinementFlags::Value)
                {
                  if (local_field_type == TensorRank::Scalar)
                    {
                      // Get the values for a scalar field
                      fe_values.get_function_values(solution_vector, values);
                    }
                  else
                    {
                      // Get the magnitude of the value for vector fields
                      std::vector<dealii::Vector<number>> vector_values(
                        num_quad_points,
                        dealii::Vector<number>(dim));
                      fe_values.get_function_values(solution_vector, vector_values);
                      for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                        {
                          values[q_point] = vector_values[q_point].l2_norm();
                        }
                    }

                  // Check if any of the quadrature points meet the refinement criterion
                  for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                    {
                      if (field_criterion.value_in_open_range(values[q_point]))
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
              if (field_criterion.criterion & RefinementFlags::Gradient)
                {
                  if (local_field_type == TensorRank::Scalar)
                    {
                      // Get the magnitude of the gradient for a scalar field
                      std::vector<dealii::Tensor<1, dim, number>> scalar_gradients(
                        num_quad_points);
                      fe_values.get_function_gradients(solution_vector, scalar_gradients);
                      for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                        {
                          values[q_point] = scalar_gradients[q_point].norm();
                        }
                    }
                  else
                    {
                      std::vector<std::vector<dealii::Tensor<1, dim, number>>>
                        vector_gradients(num_quad_points,
                                         std::vector<dealii::Tensor<1, dim, number>>(
                                           dim));
                      fe_values.get_function_gradients(solution_vector, vector_gradients);
                      for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                        {
                          dealii::Vector<number> vector_gradient_component_magnitude(dim);
                          for (unsigned int dimension = 0; dimension < dim; dimension++)
                            {
                              vector_gradient_component_magnitude[dimension] =
                                vector_gradients[q_point][dimension].norm();
                            }
                          values[q_point] = vector_gradient_component_magnitude.l2_norm();
                        }
                    }

                  // Check if any of the quadrature points meet the refinement criterion
                  for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                    {
                      if (field_criterion.gradient_magnitude_above_threshold(
                            values[q_point]))
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

          DEBUG_ASSERT(
            cell->level() > 0,
            "Cell refinement level is less than one, which will lead to underflow",
            cell->level());

          const auto cell_refinement = (unsigned int) cell->level();
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

template <unsigned int dim, unsigned int degree, typename number>
bool
RefinementManager<dim, degree, number>::mark_cells_for_refinement()
{
  Timer::Scope timer_scope("mark_cells_for_refinement()");

  bool any_cell_marked = false;
  for (const auto &cell : solve_context->get_triangulation_manager()
                            .get_triangulation()
                            .active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          const unsigned int cell_refinement = cell->level();
          if (std::any_of(
                marker_functions.begin(),
                marker_functions.end(),
                [&](const std::shared_ptr<const CellMarkerBase<dim>> &marker_function)
                {
                  return marker_function->flag(*cell,
                                               solve_context->get_simulation_timer());
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

template <unsigned int dim, unsigned int degree, typename number>
void
RefinementManager<dim, degree, number>::refine_grid(
  std::vector<std::shared_ptr<SolverBase<dim, degree, number>>> &solvers)
{
  Timer::Scope timer_scope("refine_grid()");

  TriangulationManager<dim> &triangulation_manager =
    solve_context->get_triangulation_manager();
  DoFManager<dim, degree>                &dof_manager = solve_context->get_dof_manager();
  ConstraintManager<dim, degree, number> &constraint_manager =
    solve_context->get_constraint_manager();
  MatrixFreeManager<dim, number> &matrix_free_manager =
    solve_context->get_matrix_free_manager();

  // Update ghosts of all fields.
  for (auto &solver : solvers)
    {
      solver->update_ghosts();
    }

  {
    Timer::Scope scope("prepare");
    // Prepare for grid refinement
    triangulation_manager.prepare_for_grid_refinement();
    for (auto &solver : solvers)
      {
        solver->prepare_for_solution_transfer();
      }
  }

  {
    Timer::Scope scope("execute");
    // Execute grid refinement
    triangulation_manager.execute_grid_refinement();
  }

  {
    Timer::Scope scope("reinit");
    // Redistribute DoFs and reinit the solvers
    if (triangulation_manager.has_mg())
      {
        triangulation_manager.init_mg();
      }
    dof_manager.reinit(triangulation_manager, dof_manager.has_mg());
    constraint_manager.reinit(solve_context->get_field_attributes());
    matrix_free_manager.reinit(dof_manager, constraint_manager);

    // Reinit solutions, apply constraints, then solution transfer
    for (auto &solver : solvers)
      {
        solver->reinit();
        solver->execute_solution_transfer();
      }
  }
}

#include "core/refinement_manager.inst"

PRISMS_PF_END_NAMESPACE
