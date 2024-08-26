// methods to apply boundary conditons

#include "../../include/matrixFreePDE.h"
#include "../../include/nonUniformDirichletBC.h"
#include "../../include/varBCs.h"
#include "../../include/vectorBCFunction.h"

// =================================================================================
// Methods to apply non-zero Neumann BCs
// =================================================================================
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::applyNeumannBCs()
{
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  // NOTE: Currently this function doesn't work and it's call is commented out
  // in solveIncrement. The result is off by almost exactly a factor of 100,000.
  // I don't know what the issue is.
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  // Check to the BC for the current field
  unsigned int starting_BC_list_index = 0;
  for (unsigned int i = 0; i < currentFieldIndex; i++)
    {
      if (userInputs.var_type[i] == SCALAR)
        {
          starting_BC_list_index++;
        }
      else
        {
          starting_BC_list_index += dim;
        }
    }

  if (userInputs.var_type[currentFieldIndex] == SCALAR)
    {
      for (unsigned int direction = 0; direction < 2 * dim; direction++)
        {
          if (userInputs.BC_list[starting_BC_list_index].var_BC_type[direction] ==
              NEUMANN)
            {
              typename DoFHandler<dim>::active_cell_iterator cell = dofHandlersSet[0]
                                                                      ->begin_active(),
                                                             endc =
                                                               dofHandlersSet[0]->end();
              FESystem<dim>         *fe = FESet[currentFieldIndex];
              QGaussLobatto<dim - 1> face_quadrature_formula(degree + 1);
              FEFaceValues<dim>      fe_face_values(*fe,
                                               face_quadrature_formula,
                                               update_values | update_JxW_values);
              const unsigned int     n_face_q_points = face_quadrature_formula.size(),
                                 dofs_per_cell       = fe->dofs_per_cell;
              Vector<double>                       cell_rhs(dofs_per_cell);
              std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

              // Loop over each face on a boundary
              for (; cell != endc; ++cell)
                {
                  for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
                    {
                      if (cell->face(f)->at_boundary())
                        {
                          if (cell->face(f)->boundary_id() == direction)
                            {
                              fe_face_values.reinit(cell, f);
                              cell_rhs = 0.0;
                              for (unsigned int q_point = 0; q_point < n_face_q_points;
                                   ++q_point)
                                {
                                  double neumann_value =
                                    userInputs.BC_list[starting_BC_list_index]
                                      .var_BC_val[direction];
                                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                                    {
                                      cell_rhs(i) +=
                                        (neumann_value *
                                         fe_face_values.shape_value(i, q_point) *
                                         fe_face_values.JxW(q_point));
                                    }
                                }
                              cell->get_dof_indices(local_dof_indices);
                              // assemble
                              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                                {
                                  (*(residualSet[currentFieldIndex]))
                                    [local_dof_indices[i]] += cell_rhs(i);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

// =================================================================================
// Methods to apply non-zero Dirichlet BCs
// =================================================================================
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::applyDirichletBCs()
{
  // First, get the variable index of the current field
  unsigned int starting_BC_list_index = 0;

  for (unsigned int i = 0; i < currentFieldIndex; i++)
    {
      if (userInputs.var_type[i] == SCALAR)
        {
          starting_BC_list_index++;
        }
      else
        {
          starting_BC_list_index += dim;
        }
    }

  if (userInputs.var_type[currentFieldIndex] == SCALAR)
    {
      for (unsigned int direction = 0; direction < 2 * dim; direction++)
        {
          if (userInputs.BC_list[starting_BC_list_index].var_BC_type[direction] ==
              DIRICHLET)
            {
              VectorTools::interpolate_boundary_values(
                *dofHandlersSet[currentFieldIndex],
                direction,
                Functions::ConstantFunction<dim>(
                  userInputs.BC_list[starting_BC_list_index].var_BC_val[direction],
                  1),
                *(AffineConstraints<double> *)
                  constraintsDirichletSet[currentFieldIndex]);
            }
          else if (userInputs.BC_list[starting_BC_list_index].var_BC_type[direction] ==
                   NON_UNIFORM_DIRICHLET)
            {
              VectorTools::interpolate_boundary_values(
                *dofHandlersSet[currentFieldIndex],
                direction,
                NonUniformDirichletBC<dim, degree>(currentFieldIndex,
                                                   direction,
                                                   currentTime,
                                                   this),
                *(AffineConstraints<double> *)
                  constraintsDirichletSet[currentFieldIndex]);
            }
        }
    }
  else
    {
      for (unsigned int direction = 0; direction < 2 * dim; direction++)
        {
          std::vector<double> BC_values;
          for (unsigned int component = 0; component < dim; component++)
            {
              BC_values.push_back(userInputs.BC_list[starting_BC_list_index + component]
                                    .var_BC_val[direction]);
            }

          std::vector<bool> mask;
          for (unsigned int component = 0; component < dim; component++)
            {
              if (userInputs.BC_list[starting_BC_list_index + component]
                    .var_BC_type[direction] == DIRICHLET)
                {
                  mask.push_back(true);
                }
              else
                {
                  mask.push_back(false);
                }
            }

          VectorTools::interpolate_boundary_values(
            *dofHandlersSet[currentFieldIndex],
            direction,
            vectorBCFunction<dim>(BC_values),
            *(AffineConstraints<double> *) constraintsDirichletSet[currentFieldIndex],
            mask);

          // Mask again, this time for non-uniform Dirichlet BCs
          mask.clear();
          for (unsigned int component = 0; component < dim; component++)
            {
              if (userInputs.BC_list[starting_BC_list_index + component]
                    .var_BC_type[direction] == NON_UNIFORM_DIRICHLET)
                {
                  mask.push_back(true);
                }
              else
                {
                  mask.push_back(false);
                }
            }

          // VectorTools::interpolate_boundary_values (*dofHandlersSet[currentFieldIndex],\
				//   direction, NonUniformDirichletBC<dim,degree>(currentFieldIndex,direction,currentTime,this), *(AffineConstraints<double>*) \
				//   constraintsDirichletSet[currentFieldIndex],mask);
          VectorTools::interpolate_boundary_values(
            *dofHandlersSet[currentFieldIndex],
            direction,
            NonUniformDirichletBCVector<dim, degree>(currentFieldIndex,
                                                     direction,
                                                     currentTime,
                                                     this),
            *(AffineConstraints<double> *) constraintsDirichletSet[currentFieldIndex],
            mask);
        }
    }
}

// Based on the contents of BC_list, mark faces on the triangulation as periodic
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::setPeriodicity()
{
  std::vector<GridTools::PeriodicFacePair<
    typename parallel::distributed::Triangulation<dim>::cell_iterator>>
    periodicity_vector;
  for (int i = 0; i < dim; ++i)
    {
      bool periodic_pair = false;
      for (unsigned int field_num = 0; field_num < userInputs.BC_list.size(); field_num++)
        {
          if (userInputs.BC_list[field_num].var_BC_type[2 * i] == PERIODIC)
            {
              periodic_pair = true;
            }
        }
      if (periodic_pair == true)
        {
          GridTools::collect_periodic_faces(triangulation,
                                            /*b_id1*/ 2 * i,
                                            /*b_id2*/ 2 * i + 1,
                                            /*direction*/ i,
                                            periodicity_vector);
        }
    }

  triangulation.add_periodicity(periodicity_vector);
  pcout << "periodic facepairs: " << periodicity_vector.size() << std::endl;
}

// Set constraints to enforce periodic boundary conditions
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::setPeriodicityConstraints(
  AffineConstraints<double> *constraints,
  const DoFHandler<dim>     *dof_handler) const
{
  // First, get the variable index of the current field
  unsigned int starting_BC_list_index = 0;
  for (unsigned int i = 0; i < currentFieldIndex; i++)
    {
      if (userInputs.var_type[i] == SCALAR)
        {
          starting_BC_list_index++;
        }
      else
        {
          starting_BC_list_index += dim;
        }
    }

  std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator>>
    periodicity_vector;
  for (int i = 0; i < dim; ++i)
    {
      if (userInputs.BC_list[starting_BC_list_index].var_BC_type[2 * i] == PERIODIC)
        {
          GridTools::collect_periodic_faces(*dof_handler,
                                            /*b_id1*/ 2 * i,
                                            /*b_id2*/ 2 * i + 1,
                                            /*direction*/ i,
                                            periodicity_vector);
        }
    }
#if (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR >= 4)
  DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vector, *constraints);
#else
  DoFTools::make_periodicity_constraints<DoFHandler<dim>>(periodicity_vector,
                                                          *constraints);
#endif
}

// Determine which (if any) components of the current field have rigid body
// modes (i.e no Dirichlet BCs) if the equation is elliptic
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::getComponentsWithRigidBodyModes(
  std::vector<int> &rigidBodyModeComponents) const
{
  // Rigid body modes only matter for elliptic equations
  if (userInputs.var_eq_type[currentFieldIndex] == IMPLICIT_TIME_DEPENDENT ||
      userInputs.var_eq_type[currentFieldIndex] == TIME_INDEPENDENT)
    {
      // First, get the variable index of the current field
      unsigned int starting_BC_list_index = 0;
      for (unsigned int i = 0; i < currentFieldIndex; i++)
        {
          if (userInputs.var_type[i] == SCALAR)
            {
              starting_BC_list_index++;
            }
          else
            {
              starting_BC_list_index += dim;
            }
        }

      // Get number of components of the field
      unsigned int num_components = 1;
      if (userInputs.var_type[currentFieldIndex] == VECTOR)
        {
          num_components = dim;
        }

      // Loop over each component and determine if it has a rigid body mode
      // (i.e. no Dirichlet BCs)
      for (unsigned int component = 0; component < num_components; component++)
        {
          bool rigidBodyMode = true;
          for (unsigned int direction = 0; direction < 2 * dim; direction++)
            {
              if (userInputs.BC_list[starting_BC_list_index + component]
                    .var_BC_type[direction] == DIRICHLET)
                {
                  rigidBodyMode = false;
                }
            }
          // If the component has a rigid body mode, add it to the list
          if (rigidBodyMode == true)
            {
              rigidBodyModeComponents.push_back(component);
            }
        }
    }
}

// Set constraints to pin the solution if there are no Dirichlet BCs for a
// component of a variable in an elliptic equation
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::setRigidBodyModeConstraints(
  const std::vector<int>     rigidBodyModeComponents,
  AffineConstraints<double> *constraints,
  const DoFHandler<dim>     *dof_handler) const
{
  if (rigidBodyModeComponents.size() > 0)
    {
      // Choose the point where the constraint will be placed. Must be the
      // coordinates of a vertex.
      dealii::Point<dim> target_point; // default constructor places the point at the
                                       // origin

      unsigned int vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;

      // Loop over each locally owned cell
      typename DoFHandler<dim>::active_cell_iterator cell = dof_handler->begin_active(),
                                                     endc = dof_handler->end();

      for (; cell != endc; ++cell)
        {
          if (cell->is_locally_owned())
            {
              for (unsigned int i = 0; i < vertices_per_cell; ++i)
                {
                  // Check if the vertex is the target vertex
                  if (target_point.distance(cell->vertex(i)) < 1e-2 * cell->diameter())
                    {
                      // Loop through the list of components with rigid body
                      // modes and add an inhomogeneous constraint for each
                      for (unsigned int component_num = 0;
                           component_num < rigidBodyModeComponents.size();
                           component_num++)
                        {
                          unsigned int nodeID = cell->vertex_dof_index(i, component_num);
                          // Temporarily disabling the addition of inhomogeneous
                          // constraints constraints->add_line(nodeID);
                          // constraints->set_inhomogeneity(nodeID,0.0);
                        }
                    }
                }
            }
        }
    }
}

#include "../../include/matrixFreePDE_template_instantiations.h"
