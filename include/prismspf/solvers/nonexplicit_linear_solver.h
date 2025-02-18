#ifndef nonexplicit_linear_solver_h
#define nonexplicit_linear_solver_h

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>

#include <prismspf/config.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/nonexplicit_base.h>
#include <prismspf/user_inputs/user_input_parameters.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * Forward declaration for user-implemented PDE class.
 */
template <int dim, int degree, typename number>
class customPDE;

/**
 * \brief This class handles all linear solves.
 */
template <int dim, int degree>
class nonexplicitLinearSolver : public nonexplicitBase<dim, degree>
{
public:
  using SystemMatrixType = customPDE<dim, degree, double>;

  /**
   * \brief Constructor.
   */
  nonexplicitLinearSolver(
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
  ~nonexplicitLinearSolver() = default;

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
nonexplicitLinearSolver<dim, degree>::nonexplicitLinearSolver(
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
nonexplicitLinearSolver<dim, degree>::solve()
{
  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  for (const auto &[index, variable] : this->subset_attributes)
    {
      if (this->user_inputs.linear_solve_parameters.linear_solve.at(index)
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

PRISMS_PF_END_NAMESPACE

#endif