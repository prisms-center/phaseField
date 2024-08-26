// init() method for MatrixFreePDE class

#include <deal.II/grid/grid_generator.h>

#include "../../include/matrixFreePDE.h"
#include "../../include/varBCs.h"

// populate with fields and setup matrix free system
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::init()
{
  computing_timer.enter_subsection("matrixFreePDE: initialization");

  // creating mesh

  pcout << "creating problem mesh...\n";
  // Create the coarse mesh and mark the boundaries
  makeTriangulation(triangulation);

  // Set which (if any) faces of the triangulation are periodic
  setPeriodicity();

  // If resuming from a checkpoint, load the refined triangulation, otherwise
  // refine globally per the parameters.in file
  if (userInputs.resume_from_checkpoint)
    {
      load_checkpoint_triangulation();
    }
  else
    {
      // Do the initial global refinement
      triangulation.refine_global(userInputs.refine_factor);
    }

  // Write out the size of the computational domain and the total number of
  // elements
  if (dim < 3)
    {
      pcout << "problem dimensions: " << userInputs.domain_size[0] << "x"
            << userInputs.domain_size[1] << std::endl;
    }
  else
    {
      pcout << "problem dimensions: " << userInputs.domain_size[0] << "x"
            << userInputs.domain_size[1] << "x" << userInputs.domain_size[2] << std::endl;
    }
  pcout << "number of elements: " << triangulation.n_global_active_cells() << std::endl;
  pcout << std::endl;

  // Setup system
  pcout << "initializing matrix free object\n";
  totalDOFs = 0;
  for (typename std::vector<Field<dim>>::iterator it = fields.begin(); it != fields.end();
       ++it)
    {
      currentFieldIndex = it->index;

      char buffer[100];

      // print to std::out
      std::string var_type;
      if (it->pdetype == EXPLICIT_TIME_DEPENDENT)
        {
          var_type = "EXPLICIT_TIME_DEPENDENT";
        }
      else if (it->pdetype == IMPLICIT_TIME_DEPENDENT)
        {
          var_type = "IMPLICIT_TIME_DEPENDENT";
        }
      else if (it->pdetype == TIME_INDEPENDENT)
        {
          var_type = "TIME_INDEPENDENT";
        }
      else if (it->pdetype == AUXILIARY)
        {
          var_type = "AUXILIARY";
        }

      snprintf(buffer,
               sizeof(buffer),
               "initializing finite element space P^%u for %9s:%6s field '%s'\n",
               degree,
               var_type.c_str(),
               (it->type == SCALAR ? "SCALAR" : "VECTOR"),
               it->name.c_str());
      pcout << buffer;

      // Check if any time dependent fields present (note: I should get rid of
      // parabolicFieldIndex and ellipticFieldIndex, they only work if there is
      // at max one of each)
      if (it->pdetype == EXPLICIT_TIME_DEPENDENT)
        {
          isTimeDependentBVP  = true;
          parabolicFieldIndex = it->index;
          hasExplicitEquation = true;
        }
      else if (it->pdetype == IMPLICIT_TIME_DEPENDENT)
        {
          isTimeDependentBVP     = true;
          ellipticFieldIndex     = it->index;
          hasNonExplicitEquation = true;
          std::cerr << "PRISMS-PF Error: IMPLICIT_TIME_DEPENDENT equation "
                       "types are not currently supported"
                    << std::endl;
          abort();
        }
      else if (it->pdetype == AUXILIARY)
        {
          parabolicFieldIndex    = it->index;
          ellipticFieldIndex     = it->index;
          hasNonExplicitEquation = true;
        }
      else if (it->pdetype == TIME_INDEPENDENT)
        {
          isEllipticBVP          = true;
          ellipticFieldIndex     = it->index;
          hasNonExplicitEquation = true;
        }

      // create FESystem
      FESystem<dim> *fe;

      if (it->type == SCALAR)
        {
          fe = new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(degree + 1)), 1);
        }
      else if (it->type == VECTOR)
        {
          fe = new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(degree + 1)), dim);
        }
      else
        {
          pcout << "\nmatrixFreePDE.h: unknown field type\n";
          exit(-1);
        }
      FESet.push_back(fe);

      // distribute DOFs
      DoFHandler<dim> *dof_handler;

      dof_handler = new DoFHandler<dim>(triangulation);
      dofHandlersSet.push_back(dof_handler);
      dofHandlersSet_nonconst.push_back(dof_handler);

      dof_handler->distribute_dofs(*fe);
      totalDOFs += dof_handler->n_dofs();

      // Extract locally_relevant_dofs
      IndexSet *locally_relevant_dofs;

      locally_relevant_dofs = new IndexSet;
      locally_relevant_dofsSet.push_back(locally_relevant_dofs);
      locally_relevant_dofsSet_nonconst.push_back(locally_relevant_dofs);

      locally_relevant_dofs->clear();
      DoFTools::extract_locally_relevant_dofs(*dof_handler, *locally_relevant_dofs);

      // Create constraints
      AffineConstraints<double> *constraintsDirichlet, *constraintsOther;

      constraintsDirichlet = new AffineConstraints<double>;
      constraintsDirichletSet.push_back(constraintsDirichlet);
      constraintsDirichletSet_nonconst.push_back(constraintsDirichlet);
      constraintsOther = new AffineConstraints<double>;
      constraintsOtherSet.push_back(constraintsOther);
      constraintsOtherSet_nonconst.push_back(constraintsOther);
      valuesDirichletSet.push_back(new std::map<dealii::types::global_dof_index, double>);

      constraintsDirichlet->clear();
      constraintsDirichlet->reinit(*locally_relevant_dofs);
      constraintsOther->clear();
      constraintsOther->reinit(*locally_relevant_dofs);

      // Get hanging node constraints
      DoFTools::make_hanging_node_constraints(*dof_handler, *constraintsOther);

      // Add a constraint to fix the value at the origin to zero if all BCs are
      // zero-derivative or periodic
      std::vector<int> rigidBodyModeComponents;
      // getComponentsWithRigidBodyModes(rigidBodyModeComponents);
      // setRigidBodyModeConstraints(rigidBodyModeComponents,constraintsOther,dof_handler);

      // Get constraints for periodic BCs
      setPeriodicityConstraints(constraintsOther, dof_handler);

      // Check if Dirichlet BCs are used
      for (unsigned int i = 0; i < userInputs.BC_list.size(); i++)
        {
          for (unsigned int direction = 0; direction < 2 * dim; direction++)
            {
              if (userInputs.BC_list[i].var_BC_type[direction] == DIRICHLET)
                {
                  it->hasDirichletBCs = true;
                }
              else if (userInputs.BC_list[i].var_BC_type[direction] ==
                       NON_UNIFORM_DIRICHLET)
                {
                  it->hasnonuniformDirichletBCs = true;
                }
              else if (userInputs.BC_list[i].var_BC_type[direction] == NEUMANN)
                {
                  it->hasNeumannBCs = true;
                }
            }
        }

      // Get constraints for Dirichlet BCs
      applyDirichletBCs();

      constraintsDirichlet->close();
      constraintsOther->close();

      // Store Dirichlet BC DOF's
      valuesDirichletSet[it->index]->clear();
      for (types::global_dof_index i = 0; i < dof_handler->n_dofs(); i++)
        {
          if (locally_relevant_dofs->is_element(i))
            {
              if (constraintsDirichlet->is_constrained(i))
                {
                  (*valuesDirichletSet[it->index])[i] =
                    constraintsDirichlet->get_inhomogeneity(i);
                }
            }
        }

      snprintf(buffer,
               sizeof(buffer),
               "field '%2s' DOF : %u (Constraint DOF : %u)\n",
               it->name.c_str(),
               dof_handler->n_dofs(),
               constraintsDirichlet->n_constraints());
      pcout << buffer;
    }
  pcout << "total DOF : " << totalDOFs << std::endl;

  // Setup the matrix free object
  typename MatrixFree<dim, double>::AdditionalData additional_data;
// The member "mpi_communicator" was removed in deal.II version 8.5 but is
// required before it
#if (DEAL_II_VERSION_MAJOR < 9 && DEAL_II_VERSION_MINOR < 5)
  additional_data.mpi_communicator = MPI_COMM_WORLD;
#endif
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, double>::AdditionalData::partition_partition;
  // additional_data.tasks_parallel_scheme =
  // MatrixFree<dim,double>::AdditionalData::none;
  additional_data.mapping_update_flags =
    (update_values | update_gradients | update_JxW_values | update_quadrature_points);
  QGaussLobatto<1> quadrature(degree + 1);
  matrixFreeObject.clear();
#if (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR < 4)
  matrixFreeObject.reinit(dofHandlersSet,
                          constraintsOtherSet,
                          quadrature,
                          additional_data);
#else
  matrixFreeObject.reinit(MappingFE<dim, dim>(FE_Q<dim>(QGaussLobatto<1>(degree + 1))),
                          dofHandlersSet,
                          constraintsOtherSet,
                          quadrature,
                          additional_data);
#endif
  bool dU_scalar_init = false;
  bool dU_vector_init = false;

  // Setup solution vectors
  pcout << "initializing parallel::distributed residual and solution vectors\n";
  for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
    {
      vectorType *U, *R;

      U = new vectorType;
      R = new vectorType;
      solutionSet.push_back(U);
      residualSet.push_back(R);
      matrixFreeObject.initialize_dof_vector(*R, fieldIndex);
      *R = 0;

      matrixFreeObject.initialize_dof_vector(*U, fieldIndex);
      *U = 0;

      // Initializing temporary dU vector required for implicit solves of the
      // elliptic equation.
      if (fields[fieldIndex].pdetype == TIME_INDEPENDENT ||
          fields[fieldIndex].pdetype == IMPLICIT_TIME_DEPENDENT ||
          (fields[fieldIndex].pdetype == AUXILIARY &&
           userInputs.var_nonlinear[fieldIndex]))
        {
          if (fields[fieldIndex].type == SCALAR)
            {
              if (dU_scalar_init == false)
                {
                  matrixFreeObject.initialize_dof_vector(dU_scalar, fieldIndex);
                  dU_scalar_init = true;
                }
            }
          else
            {
              if (dU_vector_init == false)
                {
                  matrixFreeObject.initialize_dof_vector(dU_vector, fieldIndex);
                  dU_vector_init = true;
                }
            }
        }
    }

  // check if time dependent BVP and compute invM
  if (isTimeDependentBVP)
    {
      computeInvM();
    }

  // Apply the initial conditions to the solution vectors
  // The initial conditions are re-applied below in the "adaptiveRefine"
  // function so that the mesh can adapt based on the initial conditions.
  if (userInputs.resume_from_checkpoint)
    {
      load_checkpoint_fields();
    }
  else
    {
      applyInitialConditions();
    }

  // Create new solution transfer sets (needed for the "refineGrid" call, might
  // be able to move this elsewhere)
  soltransSet.clear();
  for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
    {
      soltransSet.push_back(new parallel::distributed::SolutionTransfer<dim, vectorType>(
        *dofHandlersSet_nonconst[fieldIndex]));
    }

  // Ghost the solution vectors. Also apply the constraints (if any) on the
  // solution vectors
  for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
    {
      constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
      constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
      solutionSet[fieldIndex]->update_ghost_values();
    }

  // If not resuming from a checkpoint, check and perform adaptive mesh
  // refinement, which reinitializes the system with the new mesh
  if (!userInputs.resume_from_checkpoint)
    {
      adaptiveRefine(0);
    }

  // If resuming from a checkpoint, load the proper starting increment and time
  if (userInputs.resume_from_checkpoint)
    {
      load_checkpoint_time_info();
    }

  computing_timer.leave_subsection("matrixFreePDE: initialization");
}

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::makeTriangulation(
  parallel::distributed::Triangulation<dim> &tria) const
{
  if (dim == 3)
    {
      GridGenerator::subdivided_hyper_rectangle(tria,
                                                userInputs.subdivisions,
                                                Point<dim>(),
                                                Point<dim>(userInputs.domain_size[0],
                                                           userInputs.domain_size[1],
                                                           userInputs.domain_size[2]));
    }
  else if (dim == 2)
    {
      GridGenerator::subdivided_hyper_rectangle(tria,
                                                userInputs.subdivisions,
                                                Point<dim>(),
                                                Point<dim>(userInputs.domain_size[0],
                                                           userInputs.domain_size[1]));
    }
  else
    {
      GridGenerator::subdivided_hyper_rectangle(tria,
                                                userInputs.subdivisions,
                                                Point<dim>(),
                                                Point<dim>(userInputs.domain_size[0]));
    }

  // Mark boundaries for applying the boundary conditions
  markBoundaries(tria);
}

#include "../../include/matrixFreePDE_template_instantiations.h"
