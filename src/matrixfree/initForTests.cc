// initForTests() method for MatrixFreePDE class

#include <deal.II/grid/grid_generator.h>

#include "../../include/matrixFreePDE.h"

template <int dim, int degree>
void
MatrixFreePDE<dim, degree>::initForTests(std::vector<Field<dim>> fields)
{
  // creating mesh
  std::vector<unsigned int> subdivisions;
  for (unsigned int i = 0; i < dim; i++)
    {
      subdivisions.push_back(10);
    }

  if (dim == 3)
    {
      GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                subdivisions,
                                                Point<dim>(),
                                                Point<dim>(1, 1, 1));
    }
  else if (dim == 2)
    {
      GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                subdivisions,
                                                Point<dim>(),
                                                Point<dim>(1, 1));
    }
  else
    {
      GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                subdivisions,
                                                Point<dim>(),
                                                Point<dim>(1));
    }

  // setup system
  for (auto &field : fields)
    {
      // create FESystem
      FESystem<dim> *fe = nullptr;
      if (field.type == SCALAR)
        {
          fe = new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(degree + 1)), 1);
        }
      else if (field.type == VECTOR)
        {
          fe = new FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(degree + 1)), dim);
        }
      FESet.push_back(fe);

      // distribute DOFs
      DoFHandler<dim> *dof_handler = nullptr;
      dof_handler                  = new DoFHandler<dim>(triangulation);
      dofHandlersSet.push_back(dof_handler);
      dofHandlersSet_nonconst.push_back(dof_handler);
      dof_handler->distribute_dofs(*fe);

      // extract locally_relevant_dofs
      IndexSet *locally_relevant_dofs = nullptr;
      locally_relevant_dofs           = new IndexSet;
      locally_relevant_dofsSet.push_back(locally_relevant_dofs);
      locally_relevant_dofsSet_nonconst.push_back(locally_relevant_dofs);
      locally_relevant_dofs->clear();
      DoFTools::extract_locally_relevant_dofs(*dof_handler, *locally_relevant_dofs);

      // create constraints
      AffineConstraints<double> *constraintsOther = nullptr;
      constraintsOther                            = new AffineConstraints<double>;
      constraintsOtherSet.push_back(constraintsOther);
      constraintsOtherSet_nonconst.push_back(constraintsOther);
      constraintsOther->clear();
      constraintsOther->reinit(*locally_relevant_dofs);
      DoFTools::make_hanging_node_constraints(*dof_handler, *constraintsOther);
      // constraintsOther->close(); // I actually don't want to close them since
      // I'll be adding to them in test_setRigidBodyModeConstraints
    }

  // //setup the matrix free object
  typename MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.tasks_parallel_scheme =
    MatrixFree<dim, double>::AdditionalData::partition_partition;
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
  // setup problem vectors
  for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
    {
      vectorType *U = nullptr;
      vectorType *R = nullptr;
      U             = new vectorType;
      R             = new vectorType;
      solutionSet.push_back(U);
      residualSet.push_back(R);
      matrixFreeObject.initialize_dof_vector(*R, 0);
      *R = 0;
      matrixFreeObject.initialize_dof_vector(*U, 0);
      *U = 0;
    }
}
