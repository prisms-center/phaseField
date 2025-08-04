
#include <prismspf/core/matrix_free_operator.h>

#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/sequential_linear_solver.h>
#include <prismspf/solvers/sequential_solver.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
SequentialLinearSolver<dim, degree, number>::SequentialLinearSolver(
  const SolverContext<dim, degree, number> &_solver_context,
  Types::Index                              _solve_priority)
  : SequentialSolver<dim, degree, number>(_solver_context,
                                          FieldSolveType::NonexplicitLinear,
                                          _solve_priority)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialLinearSolver<dim, degree, number>::init()
{
  // Call the base class init
  this->SequentialSolver<dim, degree, number>::init();

  // If the solver is empty we can just return early.
  if (this->solver_is_empty())
    {
      return;
    }

  // Init the linear solvers
  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      this->init_linear_solver(variable);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialLinearSolver<dim, degree, number>::reinit()
{
  // Call the base class reinit
  this->SequentialSolver<dim, degree, number>::reinit();

  // If the solver is empty we can just return early.
  if (this->solver_is_empty())
    {
      return;
    }

  // Reinit the linear solvers
  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      this->reinit_linear_solver(variable);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialLinearSolver<dim, degree, number>::solve()
{
  // Call the base class solve
  this->SequentialSolver<dim, degree, number>::solve();

  // If the solver is empty we can just return early.
  if (this->solver_is_empty())
    {
      return;
    }

  // Solve each field
  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      this->solve_linear_solver(variable);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SequentialLinearSolver<dim, degree, number>::print()
{
  // Print the base class information
  this->SequentialSolver<dim, degree, number>::print();
}

#include "solvers/sequential_linear_solver.inst"

PRISMS_PF_END_NAMESPACE