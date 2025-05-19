// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/mg_level_object.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/nonexplicit_base.h>
#include <prismspf/solvers/nonexplicit_self_nonlinear_solver.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
nonexplicitSelfNonlinearSolver<dim, degree>::nonexplicitSelfNonlinearSolver(
  const UserInputParameters<dim>                         &_user_inputs,
  const MatrixfreeHandler<dim>                           &_matrix_free_handler,
  const TriangulationHandler<dim>                        &_triangulation_handler,
  const InvmHandler<dim, degree>                         &_invm_handler,
  const ConstraintHandler<dim, degree>                   &_constraint_handler,
  const DofHandler<dim>                                  &_dof_handler,
  const dealii::MappingQ1<dim>                           &_mapping,
  dealii::MGLevelObject<MatrixfreeHandler<dim, float>>   &_mg_matrix_free_handler,
  SolutionHandler<dim>                                   &_solution_handler,
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
nonexplicitSelfNonlinearSolver<dim, degree>::init()
{
  this->compute_subset_attributes(FieldSolveType::NONEXPLICIT_SELF_NONLINEAR);

  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }

  this->set_initial_condition();

  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      if (this->get_user_inputs()
            .get_linear_solve_parameters()
            .get_linear_solve_parameters(index)
            .preconditioner == PreconditionerType::GMG)
        {
          gmg_solvers.emplace(
            index,
            std::make_unique<GMGSolver<dim, degree>>(this->get_user_inputs(),
                                                     variable,
                                                     this->get_matrix_free_handler(),
                                                     this->get_constraint_handler(),
                                                     this->get_triangulation_handler(),
                                                     this->get_dof_handler(),
                                                     this->get_mg_matrix_free_handler(),
                                                     this->get_solution_handler(),
                                                     this->get_pde_operator(),
                                                     pde_operator_float,
                                                     *mg_info));
          gmg_solvers.at(index)->init();
        }
      else
        {
          identity_solvers.emplace(
            index,
            std::make_unique<identitySolver<dim, degree>>(this->get_user_inputs(),
                                                          variable,
                                                          this->get_matrix_free_handler(),
                                                          this->get_constraint_handler(),
                                                          this->get_solution_handler(),
                                                          this->get_pde_operator()));
          identity_solvers.at(index)->init();
        }
    }
}

template <unsigned int dim, unsigned int degree>
inline void
nonexplicitSelfNonlinearSolver<dim, degree>::solve()
{
  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }

  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      // Skip if the field type is IMPLICIT_TIME_DEPENDENT and the current increment is 0.
      if (variable.get_pde_type() == PDEType::IMPLICIT_TIME_DEPENDENT &&
          this->get_user_inputs().get_temporal_discretization().get_current_increment() ==
            0)
        {
          continue;
        }

      bool         is_converged = true;
      unsigned int iteration    = 0;
      const auto  &step_length  = this->get_user_inputs()
                                  .get_nonlinear_solve_parameters()
                                  .get_nonlinear_solve_parameters(index)
                                  .step_length;

      while (is_converged)
        {
          is_converged = false;

          // Perform the linear solve with the step length
          if (this->get_user_inputs()
                .get_linear_solve_parameters()
                .get_linear_solve_parameters(index)
                .preconditioner == PreconditionerType::GMG)
            {
              gmg_solvers.at(index)->solve(step_length);
            }
          else
            {
              identity_solvers.at(index)->solve(step_length);
            }

          iteration++;

          if (iteration < this->get_user_inputs()
                            .get_nonlinear_solve_parameters()
                            .get_nonlinear_solve_parameters(index)
                            .max_iterations)
            {
              is_converged = true;
            }
        }
    }
}

INSTANTIATE_BI_TEMPLATE(nonexplicitSelfNonlinearSolver)

PRISMS_PF_END_NAMESPACE
