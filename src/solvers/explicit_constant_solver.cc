// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/solvers/explicit_base.h>
#include <prismspf/solvers/explicit_constant_solver.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
ExplicitConstantSolver<dim, degree>::ExplicitConstantSolver(
  const SolverContext<dim, degree> &_solver_context)
  : ExplicitBase<dim, degree>(_solver_context)
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