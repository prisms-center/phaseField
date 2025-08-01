
#include <prismspf/solvers/concurrent_explicit_postprocess_solver.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
ConcurrentExplicitPostprocessSolver<dim, degree, number>::
  ConcurrentExplicitPostprocessSolver(
    const SolverContext<dim, degree, number> &_solver_context,
    Types::Index                              _solve_priority)
  : ConcurrentSolver<dim, degree, number>(_solver_context,
                                          FieldSolveType::ExplicitPostprocess,
                                          _solve_priority)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentExplicitPostprocessSolver<dim, degree, number>::init()
{
  // Call the base class init
  this->ConcurrentSolver<dim, degree, number>::init();

  // Do nothing
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentExplicitPostprocessSolver<dim, degree, number>::reinit()
{
  // Call the base class reinit
  this->ConcurrentSolver<dim, degree, number>::reinit();

  // Do nothing
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentExplicitPostprocessSolver<dim, degree, number>::solve()
{
  // Call the base class solve
  this->ConcurrentSolver<dim, degree, number>::solve();

  // If the solver is empty we can just return early.
  if (this->solver_is_empty())
    {
      return;
    }

  // Otherwise, solve
  this->solve_explicit_equations(
    [this](std::vector<typename SolverBase<dim, degree, number>::VectorType *>       &dst,
           const std::vector<typename SolverBase<dim, degree, number>::VectorType *> &src)
    {
      this->get_system_matrix()->compute_postprocess_explicit_update(dst, src);
    });
}

template <unsigned int dim, unsigned int degree, typename number>
void
ConcurrentExplicitPostprocessSolver<dim, degree, number>::print()
{
  // Print the base class information
  this->ConcurrentSolver<dim, degree, number>::print();
}

#include "solvers/concurrent_explicit_postprocess_solver.inst"

PRISMS_PF_END_NAMESPACE