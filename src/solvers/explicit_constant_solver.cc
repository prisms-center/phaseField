// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/explicit_base.h>
#include <prismspf/solvers/explicit_constant_solver.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
ExplicitConstantSolver<dim, degree>::ExplicitConstantSolver(
  const UserInputParameters<dim>                         &_user_inputs,
  const MatrixfreeHandler<dim>                           &_matrix_free_handler,
  const InvmHandler<dim, degree>                         &_invm_handler,
  const ConstraintHandler<dim, degree>                   &_constraint_handler,
  const DofHandler<dim>                                  &_dof_handler,
  const dealii::MappingQ1<dim>                           &_mapping,
  SolutionHandler<dim>                                   &_solution_handler,
  std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator)
  : ExplicitBase<dim, degree>(_user_inputs,
                              _matrix_free_handler,
                              _invm_handler,
                              _constraint_handler,
                              _dof_handler,
                              _mapping,
                              _solution_handler,
                              std::move(_pde_operator))
{}

template <unsigned int dim, unsigned int degree>
void
ExplicitConstantSolver<dim, degree>::init()
{
  this->compute_subset_attributes(FieldSolveType::ExplicitConstant);

  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }

  this->set_initial_condition();
}

template <unsigned int dim, unsigned int degree>
void
ExplicitConstantSolver<dim, degree>::solve()
{}

INSTANTIATE_BI_TEMPLATE(ExplicitConstantSolver)

PRISMS_PF_END_NAMESPACE