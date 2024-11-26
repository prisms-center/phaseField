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
      if (var_attributes.attributes.at(i).var_type == SCALAR)
        {
          starting_BC_list_index++;
        }
      else
        {
          starting_BC_list_index += dim;
        }
    }

  if (var_attributes.attributes.at(currentFieldIndex).var_type == SCALAR)
    {
      for (unsigned int direction = 0; direction < 2 * dim; direction++)
        {
          if (userInputs.BC_list[starting_BC_list_index].var_BC_type[direction] ==
              NEUMANN)
            {
              FESystem<dim>         *fe = FESet[currentFieldIndex];
              QGaussLobatto<dim - 1> face_quadrature_formula(degree + 1);
              FEFaceValues<dim>      fe_face_values(*fe,
                                               face_quadrature_formula,
                                               update_values | update_JxW_values);
              const unsigned int     n_face_q_points = face_quadrature_formula.size();
              const unsigned int     dofs_per_cell   = fe->dofs_per_cell;
              Vector<double>         cell_rhs(dofs_per_cell);
              std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

              // Loop over each face on a boundary
              for (const auto &cell : dofHandlersSet[0]->active_cell_iterators())
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
      if (var_attributes.attributes.at(i).var_type == SCALAR)
        {
          starting_BC_list_index++;
        }
      else
        {
          starting_BC_list_index += dim;
        }
    }

  if (var_attributes.attributes.at(currentFieldIndex).var_type == SCALAR)
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

          // VectorTools::interpolate_boundary_values
          // (*dofHandlersSet[currentFieldIndex],direction,
          // NonUniformDirichletBC<dim,degree>(currentFieldIndex,direction,currentTime,this),
          // *(AffineConstraints<double>*)constraintsDirichletSet[currentFieldIndex],mask);
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
      if (periodic_pair)
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
      if (var_attributes.attributes.at(i).var_type == SCALAR)
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

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::set_rigid_body_mode_constraints(
  AffineConstraints<double> *constraints,
  const DoFHandler<dim>     *dof_handler,
  const Point<dim>           target_point) const
{
  // Determine the number of components in the field. For a scalar field this is 1, for a
  // vector dim, etc.
  unsigned int n_components = 0;
  var_attributes.attributes.at(currentFieldIndex).var_type == VECTOR ? n_components = dim
                                                                     : n_components = 1;

  // Loop over each locally owned cell
  for (const auto &cell : dof_handler->active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
            {
              // Check if the vertex is the target vertex
              if (target_point.distance(cell->vertex(i)) < 1.0e-2 * cell->diameter())
                {
                  // Loop through the number of components and add the constraint
                  for (unsigned int component = 0; component < n_components; component++)
                    {
                      unsigned int nodeID = cell->vertex_dof_index(i, component);
                      constraints->add_line(nodeID);
                      constraints->set_inhomogeneity(nodeID, 0.0);
                    }
                }
            }
        }
    }
}
