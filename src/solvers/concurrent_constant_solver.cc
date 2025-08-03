
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>

#include <prismspf/solvers/concurrent_constant_solver.h>
#include <prismspf/solvers/concurrent_solver.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
ConcurrentConstantSolver<dim, degree, number>::ConcurrentConstantSolver(
  const SolverContext<dim, degree, number> &_solver_context,
  Types::Index                              _solve_priority)
  : ConcurrentSolver<dim, degree, number>(_solver_context,
                                          FieldSolveType::ExplicitConstant,
                                          _solve_priority)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentConstantSolver<dim, degree, number>::init()
{
  // Call the base class init
  this->ConcurrentSolver<dim, degree, number>::init();

  // Do nothing
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentConstantSolver<dim, degree, number>::reinit()
{
  // Call the base class reinit
  this->ConcurrentSolver<dim, degree, number>::reinit();

  // Do nothing
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentConstantSolver<dim, degree, number>::solve()
{
  // Call the base class solve
  this->ConcurrentSolver<dim, degree, number>::solve();

  // Do nothing
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentConstantSolver<dim, degree, number>::print()
{
  // Print the base class information
  this->ConcurrentSolver<dim, degree, number>::print();
}

#include "solvers/concurrent_constant_solver.inst"

PRISMS_PF_END_NAMESPACE