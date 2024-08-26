#ifndef INCLUDE_ADAPTIVEREFINEMENT_H_
#define INCLUDE_ADAPTIVEREFINEMENT_H_

#include <deal.II/base/quadrature.h>
#include <deal.II/base/timer.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/lac/la_parallel_vector.h>

#include "fields.h"
#include "userInputParameters.h"

#ifndef vectorType
typedef dealii::LinearAlgebra::distributed::Vector<double> vectorType;
#endif

using namespace dealii;

/**
 * This class deals with adaptive refinement. Add more comments later
 */
template <int dim, int degree>
class adaptiveRefinement
{
public:
  adaptiveRefinement(
    const userInputParameters<dim>            &_userInputs,
    parallel::distributed::Triangulation<dim> &_triangulation,
    std::vector<Field<dim>>                   &_fields,
    std::vector<vectorType *>                 &_solutionSet,
    std::vector<parallel::distributed::SolutionTransfer<dim, vectorType> *> &_soltransSet,
    std::vector<FESystem<dim> *>                                            &_FESet,
    std::vector<DoFHandler<dim> *>                 &_dofHandlersSet_nonconst,
    std::vector<const AffineConstraints<double> *> &_constraintsDirichletSet,
    std::vector<const AffineConstraints<double> *> &_constraintsOtherSet);

  // Adaptive refinement
  void
  adaptiveRefine(unsigned int _currentIncrement);

  // Method that refines the triangulation
  void
  refineGrid();

  // Current increment
  unsigned int currentIncrement;

protected:
  // Adaptive refinement criterion
  void
  adaptiveRefineCriterion();

private:
  userInputParameters<dim> userInputs;

  parallel::distributed::Triangulation<dim> &triangulation;

  std::vector<Field<dim>> &fields;

  std::vector<vectorType *> &solutionSet;

  std::vector<parallel::distributed::SolutionTransfer<dim, vectorType> *> &soltransSet;

  std::vector<FESystem<dim> *> &FESet;

  std::vector<DoFHandler<dim> *> &dofHandlersSet_nonconst;

  std::vector<const AffineConstraints<double> *> &constraintsDirichletSet;

  std::vector<const AffineConstraints<double> *> &constraintsOtherSet;
};

template <int dim, int degree>
adaptiveRefinement<dim, degree>::adaptiveRefinement(
  const userInputParameters<dim>                                          &_userInputs,
  parallel::distributed::Triangulation<dim>                               &_triangulation,
  std::vector<Field<dim>>                                                 &_fields,
  std::vector<vectorType *>                                               &_solutionSet,
  std::vector<parallel::distributed::SolutionTransfer<dim, vectorType> *> &_soltransSet,
  std::vector<FESystem<dim> *>                                            &_FESet,
  std::vector<DoFHandler<dim> *>                 &_dofHandlersSet_nonconst,
  std::vector<const AffineConstraints<double> *> &_constraintsDirichletSet,
  std::vector<const AffineConstraints<double> *> &_constraintsOtherSet)
  : userInputs(_userInputs)
  , triangulation(_triangulation)
  , fields(_fields)
  , solutionSet(_solutionSet)
  , soltransSet(_soltransSet)
  , FESet(_FESet)
  , dofHandlersSet_nonconst(_dofHandlersSet_nonconst)
  , constraintsDirichletSet(_constraintsDirichletSet)
  , constraintsOtherSet(_constraintsOtherSet)
{}

template <int dim, int degree>
void
adaptiveRefinement<dim, degree>::adaptiveRefine(unsigned int currentIncrement)
{
  if ((currentIncrement == 0))
    {
      adaptiveRefineCriterion();
      refineGrid();
    }
  else
    {
      // Apply constraints before remeshing
      for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
        {
          constraintsDirichletSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
          constraintsOtherSet[fieldIndex]->distribute(*solutionSet[fieldIndex]);
          solutionSet[fieldIndex]->update_ghost_values();
        }
      adaptiveRefineCriterion();
      refineGrid();
    }
}

template <int dim, int degree>
void
adaptiveRefinement<dim, degree>::adaptiveRefineCriterion()
{
  std::vector<std::vector<double>> valuesV;
  std::vector<std::vector<double>> gradientsV;

  QGaussLobatto<dim> quadrature(degree + 1);
  const unsigned int num_quad_points = quadrature.size();

  // Set the correct update flags
  bool need_value    = false;
  bool need_gradient = false;
  for (unsigned int field_index = 0; field_index < userInputs.refinement_criteria.size();
       field_index++)
    {
      if (userInputs.refinement_criteria[field_index].criterion_type == VALUE ||
          userInputs.refinement_criteria[field_index].criterion_type ==
            VALUE_AND_GRADIENT)
        {
          need_value = true;
        }
      else if (userInputs.refinement_criteria[field_index].criterion_type == GRADIENT ||
               userInputs.refinement_criteria[field_index].criterion_type ==
                 VALUE_AND_GRADIENT)
        {
          need_gradient = true;
        }
    }
  dealii::UpdateFlags update_flags;
  if (need_value && !need_gradient)
    {
      update_flags = update_values;
    }
  else if (!need_value && need_gradient)
    {
      update_flags = update_gradients;
    }
  else
    {
      update_flags = update_values | update_gradients;
    }

  FEValues<dim> fe_values(*FESet[userInputs.refinement_criteria[0].variable_index],
                          quadrature,
                          update_flags);

  std::vector<double>                         values(num_quad_points);
  std::vector<double>                         gradient_magnitudes(num_quad_points);
  std::vector<dealii::Tensor<1, dim, double>> gradients(num_quad_points);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dofHandlersSet_nonconst[userInputs.refinement_criteria[0].variable_index]
             ->begin_active(),
    endc =
      dofHandlersSet_nonconst[userInputs.refinement_criteria[0].variable_index]->end();

  typename parallel::distributed::Triangulation<dim>::active_cell_iterator t_cell =
    triangulation.begin_active();

  for (; cell != endc; ++cell)
    {
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);

          for (unsigned int field_index = 0;
               field_index < userInputs.refinement_criteria.size();
               field_index++)
            {
              if (need_value)
                {
                  fe_values.get_function_values(
                    *solutionSet[userInputs.refinement_criteria[field_index]
                                   .variable_index],
                    values);
                  valuesV.push_back(values);
                }
              if (need_gradient)
                {
                  fe_values.get_function_gradients(
                    *solutionSet[userInputs.refinement_criteria[field_index]
                                   .variable_index],
                    gradients);

                  for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                    {
                      gradient_magnitudes.at(q_point) = gradients.at(q_point).norm();
                    }

                  gradientsV.push_back(gradient_magnitudes);
                }
            }

          bool mark_refine = false;

          for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
            {
              for (unsigned int field_index = 0;
                   field_index < userInputs.refinement_criteria.size();
                   field_index++)
                {
                  if (userInputs.refinement_criteria[field_index].criterion_type ==
                        VALUE ||
                      userInputs.refinement_criteria[field_index].criterion_type ==
                        VALUE_AND_GRADIENT)
                    {
                      if ((valuesV[field_index][q_point] >
                           userInputs.refinement_criteria[field_index]
                             .value_lower_bound) &&
                          (valuesV[field_index][q_point] <
                           userInputs.refinement_criteria[field_index].value_upper_bound))
                        {
                          mark_refine = true;
                          break;
                        }
                    }
                  if (userInputs.refinement_criteria[field_index].criterion_type ==
                        GRADIENT ||
                      userInputs.refinement_criteria[field_index].criterion_type ==
                        VALUE_AND_GRADIENT)
                    {
                      if (gradientsV[field_index][q_point] >
                          userInputs.refinement_criteria[field_index]
                            .gradient_lower_bound)
                        {
                          mark_refine = true;
                          break;
                        }
                    }
                }
            }

          valuesV.clear();
          gradientsV.clear();

          // limit the maximal and minimal refinement depth of the mesh
          unsigned int current_level = t_cell->level();

          if ((mark_refine && current_level < userInputs.max_refinement_level))
            {
              cell->set_refine_flag();
            }
          else if (!mark_refine && current_level > userInputs.min_refinement_level)
            {
              cell->set_coarsen_flag();
            }
        }
      ++t_cell;
    }
}

template <int dim, int degree>
void
adaptiveRefinement<dim, degree>::refineGrid()
{
  // prepare and refine
  triangulation.prepare_coarsening_and_refinement();
  for (unsigned int fieldIndex = 0; fieldIndex < fields.size(); fieldIndex++)
    {
      soltransSet[fieldIndex]->prepare_for_coarsening_and_refinement(
        *solutionSet[fieldIndex]);
    }
  triangulation.execute_coarsening_and_refinement();
}

#endif
