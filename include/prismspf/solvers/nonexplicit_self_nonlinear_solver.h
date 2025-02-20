// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef nonexplicit_self_nonlinear_solver_h
#define nonexplicit_self_nonlinear_solver_h

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * Forward declaration for user-implemented PDE class.
 */
template <int dim, int degree, typename number>
class customPDE;

/**
 * \brief This class handles the self-nonlinear solves of a single nonexplicit field
 */
template <int dim, int degree>
class nonexplicitSelfNonlinearSolver : public nonexplicitBase<dim, degree>
{
public:
  using SystemMatrixType = customPDE<dim, degree, double>;

  /**
   * \brief Constructor.
   */
  nonexplicitSelfNonlinearSolver(
    const userInputParameters<dim>                       &_user_inputs,
    const matrixfreeHandler<dim>                         &_matrix_free_handler,
    const triangulationHandler<dim>                      &_triangulation_handler,
    const invmHandler<dim, degree>                       &_invm_handler,
    const constraintHandler<dim>                         &_constraint_handler,
    const prisms::dofHandler<dim>                        &_dof_handler,
    const dealii::MappingQ1<dim>                         &_mapping,
    dealii::MGLevelObject<matrixfreeHandler<dim, float>> &_mg_matrix_free_handler,
    solutionHandler<dim>                                 &_solution_handler);

  /**
   * \brief Destructor.
   */
  ~nonexplicitSelfNonlinearSolver() = default;

  /**
   * \brief Initialize system.
   */
  void
  init() override;

  /**
   * \brief Solve a single update step.
   */
  void
  solve() override;

private:
  /**
   * \brief Map of identity linear solvers
   */
  std::map<unsigned int, std::unique_ptr<identitySolver<dim, degree>>> identity_solvers;

  /**
   * \brief Map of geometric multigrid linear solvers
   */
  std::map<unsigned int, std::unique_ptr<GMGSolver<dim, degree>>> gmg_solvers;
};

template <int dim, int degree>
nonexplicitSelfNonlinearSolver<dim, degree>::nonexplicitSelfNonlinearSolver(
  const userInputParameters<dim>                       &_user_inputs,
  const matrixfreeHandler<dim>                         &_matrix_free_handler,
  const triangulationHandler<dim>                      &_triangulation_handler,
  const invmHandler<dim, degree>                       &_invm_handler,
  const constraintHandler<dim>                         &_constraint_handler,
  const prisms::dofHandler<dim>                        &_dof_handler,
  const dealii::MappingQ1<dim>                         &_mapping,
  dealii::MGLevelObject<matrixfreeHandler<dim, float>> &_mg_matrix_free_handler,
  solutionHandler<dim>                                 &_solution_handler)
  : nonexplicitBase<dim, degree>(_user_inputs,
                                 _matrix_free_handler,
                                 _triangulation_handler,
                                 _invm_handler,
                                 _constraint_handler,
                                 _dof_handler,
                                 _mapping,
                                 _mg_matrix_free_handler,
                                 _solution_handler)
{}

template <int dim, int degree>
inline void
nonexplicitSelfNonlinearSolver<dim, degree>::init()
{
  this->compute_subset_attributes(fieldSolveType::NONEXPLICIT_SELF_NONLINEAR);

  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  this->set_initial_condition();

  for (const auto &[index, variable] : this->subset_attributes)
    {
      if (this->user_inputs.linear_solve_parameters.linear_solve.at(index)
            .preconditioner == preconditionerType::GMG)
        {
          gmg_solvers.emplace(
            index,
            std::make_unique<GMGSolver<dim, degree>>(this->user_inputs,
                                                     variable,
                                                     this->matrix_free_handler,
                                                     this->constraint_handler,
                                                     this->triangulation_handler,
                                                     this->dof_handler,
                                                     this->mg_matrix_free_handler,
                                                     this->solution_handler));
          gmg_solvers.at(index)->init();
        }
      else
        {
          identity_solvers.emplace(
            index,
            std::make_unique<identitySolver<dim, degree>>(this->user_inputs,
                                                          variable,
                                                          this->matrix_free_handler,
                                                          this->constraint_handler,
                                                          this->solution_handler));
          identity_solvers.at(index)->init();
        }
    }
}

template <int dim, int degree>
inline void
nonexplicitSelfNonlinearSolver<dim, degree>::solve()
{
  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  for (const auto &[index, variable] : this->subset_attributes)
    {
      bool         is_converged = true;
      unsigned int iteration    = 0;
      auto        &step_length =
        this->user_inputs.nonlinear_solve_parameters.nonlinear_solve.at(index)
          .step_length;

      while (is_converged)
        {
          is_converged = false;

          // Perform the linear solve with the step length
          if (this->user_inputs.linear_solve_parameters.linear_solve.at(index)
                .preconditioner == preconditionerType::GMG)
            {
              gmg_solvers.at(index)->solve(step_length);
            }
          else
            {
              identity_solvers.at(index)->solve(step_length);
            }

          iteration++;

          if (iteration <
              this->user_inputs.nonlinear_solve_parameters.nonlinear_solve.at(index)
                .max_iterations)
            {
              is_converged = true;
            }
        }
    }
}

PRISMS_PF_END_NAMESPACE

#endif