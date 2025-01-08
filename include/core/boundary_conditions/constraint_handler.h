#ifndef constraint_handler_h
#define constraint_handler_h

#include <deal.II/base/point.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/numerics/vector_tools_boundary.h>

#include <cmath>
#include <core/boundary_conditions/nonuniform_dirichlet.h>
#include <core/exceptions.h>
#include <core/type_enums.h>
#include <core/user_inputs/user_input_parameters.h>
#include <memory>

/**
 * \brief The class handles the generation and application of boundary conditions based on
 * the user-inputs.
 */
template <int dim, fieldType field_type = fieldType::SCALAR>
class constraintHandler
{
public:
  /**
   * \brief Constructor.
   */
  constraintHandler(const boundaryParameters<dim> &_boundary_parameters,
                    const uint                    &_global_variable_index);

  /**
   * \brief Getter function for constraints (constant reference).
   */
  [[nodiscard]] dealii::AffineConstraints<double> &
  get_constraints();

  /**
   * \brief Make constraints based on the inputs of the construct.
   */
  void
  make_constraints(const dealii::MappingQ1<dim>  &mapping,
                   const dealii::DoFHandler<dim> &dof_handler);

private:
  /**
   * \brief Set the uniform dirichlet constraint.
   */
  void
  set_dirichlet_constraints();

  /**
   * \brief Set the dirichlet constraint for the pinned point.
   */
  void
  set_pinned_point(const std::pair<double, dealii::Point<dim>> &value_point_pair,
                   const dealii::DoFHandler<dim>               &dof_handler);

  std::unique_ptr<dealii::AffineConstraints<double>> constraints;

  boundaryParameters<dim> boundary_parameters;

  uint global_variable_index;
};

template <int dim, fieldType field_type>
constraintHandler<dim, field_type>::constraintHandler(
  const boundaryParameters<dim> &_boundary_parameters,
  const uint                    &_global_variable_index)
  : constraints(std::make_unique<dealii::AffineConstraints<double>>())
  , boundary_parameters(_boundary_parameters)
  , global_variable_index(_global_variable_index)
{
  constraints->clear();
}

template <int dim, fieldType field_type>
dealii::AffineConstraints<double> &
constraintHandler<dim, field_type>::get_constraints()
{
  return *constraints;
}

template <int dim, fieldType field_type>
void
constraintHandler<dim, field_type>::make_constraints(
  const dealii::MappingQ1<dim>  &mapping,
  const dealii::DoFHandler<dim> &dof_handler)
{
  // Reinitialize constraints
  constraints->reinit(dof_handler.locally_owned_dofs(),
                      dealii::DoFTools::extract_locally_relevant_dofs(dof_handler));

  // Make hanging node constraints
  dealii::DoFTools::make_hanging_node_constraints(dof_handler, *constraints);

  // First check the normal boundary conditions
  const auto &boundary_condition =
    boundary_parameters.boundary_condition_list.at(global_variable_index);

  for (const auto &[component, condition] : boundary_condition)
    {
      for (const auto &[boundary_id, boundary_type] : condition.boundary_condition_map)
        {
          if (boundary_type == boundaryType::NATURAL)
            {
              // Do nothing because they are naturally enforced.
              continue;
            }
          else if (boundary_type == boundaryType::DIRICHLET)
            {
              dealii::VectorTools::interpolate_boundary_values(
                mapping,
                dof_handler,
                boundary_id,
                dealii::Functions::ConstantFunction<dim>(condition.dirichlet_value_map.at(
                                                           boundary_id),
                                                         1),
                *constraints);
            }
          else if (boundary_type == boundaryType::PERIODIC)
            {
              // Skip boundary ids that are odd since those map to the even faces
              if (boundary_id % 2 != 0)
                {
                  continue;
                }
              // Create a vector of matched pairs that we fill and enforce upon the
              // constaints
              std::vector<dealii::GridTools::PeriodicFacePair<
                typename dealii::DoFHandler<dim>::cell_iterator>>
                periodicity_vector;

              // Determine the direction
              const uint direction = std::floor(boundary_id / dim);

              // Collect the matched pairs on the coarsest level of the mesh
              dealii::GridTools::collect_periodic_faces(dof_handler,
                                                        boundary_id,
                                                        boundary_id + 1,
                                                        direction,
                                                        periodicity_vector);

              // Set constraints
              dealii::DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vector,
                                                                       *constraints);
            }
          else if (boundary_type == boundaryType::NEUMANN)
            {
              Assert(false, FeatureNotImplemented("Neumann boundary conditions"));
            }
          else if (boundary_type == boundaryType::NON_UNIFORM_DIRICHLET)
            {
              dealii::VectorTools::interpolate_boundary_values(
                mapping,
                dof_handler,
                boundary_id,
                nonuniformDirichlet<dim, field_type>(global_variable_index, boundary_id),
                *constraints);
            }
          else if (boundary_type == boundaryType::NON_UNIFORM_NEUMANN)
            {
              Assert(false,
                     FeatureNotImplemented("Nonuniform neumann boundary conditions"));
            }
        }
    }

  // Second check for pinned points, if they exist
  if (boundary_parameters.pinned_point_list.find(global_variable_index) !=
      boundary_parameters.pinned_point_list.end())
    {
      set_pinned_point(boundary_parameters.pinned_point_list.at(global_variable_index),
                       dof_handler);
    }

  // Close constraints
  constraints->close();
}

template <int dim, fieldType field_type>
void
constraintHandler<dim, field_type>::set_pinned_point(
  const std::pair<double, dealii::Point<dim>> &value_point_pair,
  const dealii::DoFHandler<dim>               &dof_handler)
{
  Assert(field_type == fieldType::VECTOR,
         FeatureNotImplemented("Pinned points for vector fields"));

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int i = 0; i < dealii::GeometryInfo<dim>::vertices_per_cell; ++i)
            {
              // Check if the vertex is the target vertex
              if (value_point_pair.second.distance(cell->vertex(i)) <
                  1.0e-2 * cell->diameter())
                {
                  unsigned int nodeID = cell->vertex_dof_index(i, 0);
                  constraints->add_line(nodeID);
                  constraints->set_inhomogeneity(nodeID, value_point_pair.first);
                }
            }
        }
    }
}

#endif