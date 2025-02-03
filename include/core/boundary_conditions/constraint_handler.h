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
#include <vector>

/**
 * \brief The class handles the generation and application of boundary conditions based on
 * the user-inputs.
 */
template <int dim>
class constraintHandler
{
public:
  /**
   * \brief Constructor.
   */
  constraintHandler(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Getter function for the constraints.
   */
  [[nodiscard]] std::vector<const dealii::AffineConstraints<double> *>
  get_constraints();

  /**
   * \brief Make constraints based on the inputs of the construct.
   */
  void
  make_constraints(const dealii::Mapping<dim>                   &mapping,
                   const std::vector<dealii::DoFHandler<dim> *> &dof_handlers);

private:
  /**
   * \brief Make the constrainst for a single index.
   */
  void
  make_constraint(const dealii::Mapping<dim>    &mapping,
                  const dealii::DoFHandler<dim> *dof_handler,
                  const uint                    &index);

  /**
   * \brief Set the dirichlet constraint for the pinned point.
   */
  void
  set_pinned_point(const dealii::DoFHandler<dim> *dof_handler, const uint &index);

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> &user_inputs;

  /**
   * \brief Constraints.
   */
  std::vector<dealii::AffineConstraints<double>> constraints;
};

template <int dim>
constraintHandler<dim>::constraintHandler(const userInputParameters<dim> &_user_inputs)
  : user_inputs(_user_inputs)
  , constraints(user_inputs.var_attributes.size(), dealii::AffineConstraints<double>())
{}

template <int dim>
std::vector<const dealii::AffineConstraints<double> *>
constraintHandler<dim>::get_constraints()
{
  std::vector<const dealii::AffineConstraints<double> *> temp;
  temp.reserve(constraints.size());
  for (const auto &constraint : constraints)
    {
      temp.push_back(&constraint);
    }
  return temp;
}

template <int dim>
inline void
constraintHandler<dim>::make_constraints(
  const dealii::Mapping<dim>                   &mapping,
  const std::vector<dealii::DoFHandler<dim> *> &dof_handlers)
{
  for (const auto &[index, variable] : user_inputs.var_attributes)
    {
      make_constraint(mapping, dof_handlers.at(index), index);
    }
}

template <int dim>
inline void
constraintHandler<dim>::make_constraint(const dealii::Mapping<dim>    &mapping,
                                        const dealii::DoFHandler<dim> *dof_handler,
                                        const uint                    &index)
{
  // Clear constraints
  constraints.at(index).clear();

  // Reinitialize constraints
  constraints.at(index).reinit(dof_handler->locally_owned_dofs(),
                               dealii::DoFTools::extract_locally_relevant_dofs(
                                 *dof_handler));

  // Make hanging node constraints
  dealii::DoFTools::make_hanging_node_constraints(*dof_handler, constraints.at(index));

  // First check the normal boundary conditions
  const auto &boundary_condition =
    user_inputs.boundary_parameters.boundary_condition_list.at(index);

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
              Assert(user_inputs.var_attributes.at(index).field_type == fieldType::VECTOR,
                     FeatureNotImplemented("Dirichlet vector fields"));

              dealii::VectorTools::interpolate_boundary_values(
                mapping,
                *dof_handler,
                boundary_id,
                dealii::Functions::ConstantFunction<dim>(condition.dirichlet_value_map.at(
                                                           boundary_id),
                                                         1),
                constraints.at(index));
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
              dealii::GridTools::collect_periodic_faces(*dof_handler,
                                                        boundary_id,
                                                        boundary_id + 1,
                                                        direction,
                                                        periodicity_vector);

              // Set constraints
              dealii::DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vector,
                                                                       constraints.at(
                                                                         index));
            }
          else if (boundary_type == boundaryType::NEUMANN)
            {
              Assert(false, FeatureNotImplemented("Neumann boundary conditions"));
            }
          else if (boundary_type == boundaryType::NON_UNIFORM_DIRICHLET)
            {
              Assert(user_inputs.var_attributes.at(index).field_type == fieldType::VECTOR,
                     FeatureNotImplemented("Nonuniform dirichlet vector fields"));

              dealii::VectorTools::interpolate_boundary_values(
                mapping,
                *dof_handler,
                boundary_id,
                nonuniformDirichlet<dim, fieldType::SCALAR>(index, boundary_id),
                constraints.at(index));
            }
          else if (boundary_type == boundaryType::NON_UNIFORM_NEUMANN)
            {
              Assert(false,
                     FeatureNotImplemented("Nonuniform neumann boundary conditions"));
            }
        }
    }

  // Second check for pinned points, if they exist
  if (user_inputs.boundary_parameters.pinned_point_list.find(index) !=
      user_inputs.boundary_parameters.pinned_point_list.end())
    {
      set_pinned_point(dof_handler, index);
    }

  // Close constraints
  constraints.at(index).close();
}

template <int dim>
inline void
constraintHandler<dim>::set_pinned_point(const dealii::DoFHandler<dim> *dof_handler,
                                         const uint                    &index)
{
  Assert(user_inputs.var_attributes.at(index).field_type == fieldType::VECTOR,
         FeatureNotImplemented("Pinned points for vector fields"));

  const auto &value_point_pair =
    user_inputs.boundary_parameters.pinned_point_list.at(index);

  for (const auto &cell : dof_handler->active_cell_iterators())
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
                  constraints.at(index).add_line(nodeID);
                  constraints.at(index).set_inhomogeneity(nodeID, value_point_pair.first);
                }
            }
        }
    }
}

#endif