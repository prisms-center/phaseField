// reinit() method for MatrixFreePDE class

#include "../../include/matrixFreePDE.h"

// populate with fields and setup matrix free system
template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::reinit()
{
  computing_timer.enter_subsection("matrixFreePDE: reinitialization");

  // setup system
  pcout << "Reinitializing matrix free object\n";
  totalDOFs = 0;
  for (typename std::vector<Field<dim>>::iterator it = fields.begin(); it != fields.end();
       ++it)
    {
      currentFieldIndex = it->index;

      char buffer[100];

      // create FESystem
      FESystem<dim> *fe;
      fe = FESet.at(it->index);

      // distribute DOFs
      DoFHandler<dim> *dof_handler;
      dof_handler = dofHandlersSet_nonconst.at(it->index);

      dof_handler->distribute_dofs(*fe);
      totalDOFs += dof_handler->n_dofs();

      // extract locally_relevant_dofs
      IndexSet *locally_relevant_dofs;
      locally_relevant_dofs = locally_relevant_dofsSet_nonconst.at(it->index);

      locally_relevant_dofs->clear();
      DoFTools::extract_locally_relevant_dofs(*dof_handler, *locally_relevant_dofs);

      // create constraints
      AffineConstraints<double> *constraintsDirichlet, *constraintsOther;

      constraintsDirichlet = constraintsDirichletSet_nonconst.at(it->index);
      constraintsOther     = constraintsOtherSet_nonconst.at(it->index);

      constraintsDirichlet->clear();
      constraintsDirichlet->reinit(*locally_relevant_dofs);
      constraintsOther->clear();
      constraintsOther->reinit(*locally_relevant_dofs);

      // Get hanging node constraints
      DoFTools::make_hanging_node_constraints(*dof_handler, *constraintsOther);

      // Add a constraint to fix the value at the origin to zero if all BCs are
      // zero-derivative or periodic
      std::vector<int> rigidBodyModeComponents;
      getComponentsWithRigidBodyModes(rigidBodyModeComponents);
      setRigidBodyModeConstraints(rigidBodyModeComponents, constraintsOther, dof_handler);

      // Get constraints for periodic BCs
      setPeriodicityConstraints(constraintsOther, dof_handler);

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
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, double>::AdditionalData::partition_partition;
  // additional_data.tasks_parallel_scheme =
  // MatrixFree<dim,double>::AdditionalData::none;
  // additional_data.tasks_block_size = 1; // This improves performance for
  // small runs, not sure about larger runs
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
      vectorType *U;

      U = solutionSet.at(fieldIndex);

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

  // Compute invM in PDE is a time-dependent BVP
  if (isTimeDependentBVP)
    {
      computeInvM();
    }

  // Transfer solution from previous mesh
  for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
    {
      // interpolate and clear used solution transfer sets
      soltransSet[fieldIndex]->interpolate(*solutionSet[fieldIndex]);
      delete soltransSet[fieldIndex];

      // reset residual vector
      vectorType *R = residualSet.at(fieldIndex);
      matrixFreeObject.initialize_dof_vector(*R, fieldIndex);
      *R = 0;
    }

  // Create new solution transfer sets
  soltransSet.clear();
  for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
    {
      soltransSet.push_back(new parallel::distributed::SolutionTransfer<dim, vectorType>(
        *dofHandlersSet_nonconst[fieldIndex]));
    }

  // If remeshing at the zeroth time step, re-apply initial conditions so the
  // starting values are correct on the refined mesh
  if (currentIncrement == 0 && !userInputs.load_grain_structure)
    {
      applyInitialConditions();
    }

  // Ghost the solution vectors. Also apply the Dirichet BC's (if any) on the
  // solution vectors
  for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
    {
      constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
      constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
      solutionSet[fieldIndex]->update_ghost_values();
    }

  computing_timer.leave_subsection("matrixFreePDE: reinitialization");
}

#include "../../include/matrixFreePDE_template_instantiations.h"
