// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef explicit_constant_solver_h
#define explicit_constant_solver_h

#include <prismspf/config.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/solvers/explicit_base.h>
#include <prismspf/user_inputs/user_input_parameters.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class handles the explicit solves of all constant fields
 */
template <int dim, int degree>
class explicitConstantSolver : public explicitBase<dim, degree>
{
public:
  /**
   * \brief Constructor.
   */
  explicitConstantSolver(const userInputParameters<dim> &_user_inputs,
                         const matrixfreeHandler<dim>   &_matrix_free_handler,
                         const invmHandler<dim, degree> &_invm_handler,
                         const constraintHandler<dim>   &_constraint_handler,
                         const prisms::dofHandler<dim>  &_dof_handler,
                         const dealii::MappingQ1<dim>   &_mapping,
                         solutionHandler<dim>           &_solution_handler);

  /**
   * \brief Destructor.
   */
  ~explicitConstantSolver() = default;

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
};

template <int dim, int degree>
explicitConstantSolver<dim, degree>::explicitConstantSolver(
  const userInputParameters<dim> &_user_inputs,
  const matrixfreeHandler<dim>   &_matrix_free_handler,
  const invmHandler<dim, degree> &_invm_handler,
  const constraintHandler<dim>   &_constraint_handler,
  const prisms::dofHandler<dim>  &_dof_handler,
  const dealii::MappingQ1<dim>   &_mapping,
  solutionHandler<dim>           &_solution_handler)
  : explicitBase<dim, degree>(_user_inputs,
                              _matrix_free_handler,
                              _invm_handler,
                              _constraint_handler,
                              _dof_handler,
                              _mapping,
                              _solution_handler)
{}

template <int dim, int degree>
inline void
explicitConstantSolver<dim, degree>::init()
{
  this->compute_subset_attributes(fieldSolveType::EXPLICIT_CONSTANT);

  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  this->set_initial_condition();
}

template <int dim, int degree>
inline void
explicitConstantSolver<dim, degree>::solve()
{}

PRISMS_PF_END_NAMESPACE

#endif