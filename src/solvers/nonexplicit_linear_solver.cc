// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/mg_level_object.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/nonexplicit_base.h>
#include <prismspf/solvers/nonexplicit_linear_solver.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
nonexplicitLinearSolver<dim, degree>::nonexplicitLinearSolver(
  const userInputParameters<dim>                         &_user_inputs,
  const matrixfreeHandler<dim>                           &_matrix_free_handler,
  const triangulationHandler<dim>                        &_triangulation_handler,
  const invmHandler<dim, degree>                         &_invm_handler,
  const constraintHandler<dim>                           &_constraint_handler,
  const dofHandler<dim>                                  &_dof_handler,
  const dealii::MappingQ1<dim>                           &_mapping,
  dealii::MGLevelObject<matrixfreeHandler<dim, float>>   &_mg_matrix_free_handler,
  solutionHandler<dim>                                   &_solution_handler,
  std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator,
  std::shared_ptr<const PDEOperator<dim, degree, float>>  _pde_operator_float,
  const MGInfo<dim>                                      &_mg_info)
  : nonexplicitBase<dim, degree>(_user_inputs,
                                 _matrix_free_handler,
                                 _triangulation_handler,
                                 _invm_handler,
                                 _constraint_handler,
                                 _dof_handler,
                                 _mapping,
                                 _mg_matrix_free_handler,
                                 _solution_handler,
                                 std::move(_pde_operator))
  , pde_operator_float(std::move(_pde_operator_float))
  , mg_info(&_mg_info)
{}

template <unsigned int dim, unsigned int degree>
inline void
nonexplicitLinearSolver<dim, degree>::init()
{
  this->compute_subset_attributes(fieldSolveType::NONEXPLICIT_LINEAR);

  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  this->set_initial_condition();

  for (const auto &[index, variable] : this->subset_attributes)
    {
      if (this->user_inputs->linear_solve_parameters.linear_solve.at(index)
            .preconditioner == preconditionerType::GMG)
        {
          gmg_solvers.emplace(
            index,
            std::make_unique<GMGSolver<dim, degree>>(*this->user_inputs,
                                                     variable,
                                                     *this->matrix_free_handler,
                                                     *this->constraint_handler,
                                                     *this->triangulation_handler,
                                                     *this->dof_handler,
                                                     *this->mg_matrix_free_handler,
                                                     *this->solution_handler,
                                                     this->pde_operator,
                                                     pde_operator_float,
                                                     *mg_info));
          gmg_solvers.at(index)->init();
        }
      else
        {
          identity_solvers.emplace(
            index,
            std::make_unique<identitySolver<dim, degree>>(*this->user_inputs,
                                                          variable,
                                                          *this->matrix_free_handler,
                                                          *this->constraint_handler,
                                                          *this->solution_handler,
                                                          this->pde_operator));
          identity_solvers.at(index)->init();
        }
    }
}

template <unsigned int dim, unsigned int degree>
inline void
nonexplicitLinearSolver<dim, degree>::solve()
{
  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  for (const auto &[index, variable] : this->subset_attributes)
    {
      // Skip if the field type is IMPLICIT_TIME_DEPENDENT and the current increment is 0.
      if (variable.pde_type == PDEType::IMPLICIT_TIME_DEPENDENT &&
          this->user_inputs->temporal_discretization.get_current_increment() == 0)
        {
          continue;
        }

      if (this->user_inputs->linear_solve_parameters.linear_solve.at(index)
            .preconditioner == preconditionerType::GMG)
        {
          gmg_solvers.at(index)->solve();
        }
      else
        {
          identity_solvers.at(index)->solve();
        }
    }
}

INSTANTIATE_BI_TEMPLATE(nonexplicitLinearSolver)

PRISMS_PF_END_NAMESPACE
