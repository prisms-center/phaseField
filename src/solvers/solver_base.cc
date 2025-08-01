
#include <prismspf/solvers/solver_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
SolverBase<dim, degree, number>::SolverBase(
  const SolverContext<dim, degree> &_solver_context,
  const FieldSolveType             &_field_solve_type,
  Types::Index                      _solve_priority)
  : solver_context(std::make_shared<SolverContext<dim, degree>>(_solver_context))
  , field_solve_type(_field_solve_type)
  , solve_priority(_solve_priority)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverBase<dim, degree, number>::init()
{
  // Update the subset of variable attributes
  update_subset_attributes(field_solve_type, solve_priority);

  // If the subset attribute is empty return early
  if (solver_is_empty())
    {
      ConditionalOStreams::pout_base() << "  no fields for this solver exist\n"
                                       << std::flush;
      return;
    }

  // Set the initial condition
  set_initial_condition();

  // Apply constraints. This part is neccessary so they are taken into account for
  // adaptive meshing
  for (const auto &[index, variable] : subset_attributes)
    {
      get_constraint_handler().get_constraint(index).distribute(
        *(get_solution_handler().get_solution_vector(index, DependencyType::Normal)));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverBase<dim, degree, number>::reinit()
{
  // If the subset attribute is empty return early
  if (solver_is_empty())
    {
      return;
    }

  // Apply constraints. This part is neccessary so they are taken into account for
  // adaptive meshing
  for (const auto &[index, variable] : subset_attributes)
    {
      get_constraint_handler().get_constraint(index).distribute(
        *(get_solution_handler().get_solution_vector(index, DependencyType::Normal)));
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverBase<dim, degree, number>::solve()
{
  // If the subset attribute is empty return early
  if (solver_is_empty())
    {
      return;
    }

  // Do nothing
}

template <unsigned int dim, unsigned int degree, typename number>
void
SolverBase<dim, degree, number>::print()
{}

#include "solvers/solver_base.inst"

PRISMS_PF_END_NAMESPACE