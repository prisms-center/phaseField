// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>

#include <prismspf/solvers/linear_solver_base.h>

#include <prismspf/config.h>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Class that handles the assembly and solving of a field with the identity
 * preconditioner (no preconditioner)
 */
template <int dim, int degree>
class identitySolver : public linearSolverBase<dim, degree>
{
public:
  using SystemMatrixType = matrixFreeOperator<dim, degree, double>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<double>;

  /**
   * \brief Constructor.
   */
  identitySolver(const userInputParameters<dim> &_user_inputs,
                 const variableAttributes       &_variable_attributes,
                 const matrixfreeHandler<dim>   &_matrix_free_handler,
                 const constraintHandler<dim>   &_constraint_handler,
                 solutionHandler<dim>           &_solution_handler,
                 std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator);

  /**
   * \brief Destructor.
   */
  ~identitySolver() override = default;

  /**
   * \brief Initialize the system.
   */
  void
  init() override;

  /**
   * \brief Reinitialize the system.
   */
  void
  reinit() override;

  /**
   * \brief Solve the system Ax=b.
   */
  void
  solve(const double &step_length = 1.0) override;
};

template <int dim, int degree>
identitySolver<dim, degree>::identitySolver(
  const userInputParameters<dim>                         &_user_inputs,
  const variableAttributes                               &_variable_attributes,
  const matrixfreeHandler<dim>                           &_matrix_free_handler,
  const constraintHandler<dim>                           &_constraint_handler,
  solutionHandler<dim>                                   &_solution_handler,
  std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator)
  : linearSolverBase<dim, degree>(_user_inputs,
                                  _variable_attributes,
                                  _matrix_free_handler,
                                  _constraint_handler,
                                  _solution_handler,
                                  _pde_operator)
{}

template <int dim, int degree>
inline void
identitySolver<dim, degree>::init()
{
  this->system_matrix->clear();
  this->system_matrix->initialize(this->matrix_free_handler->get_matrix_free());
  this->update_system_matrix->clear();
  this->update_system_matrix->initialize(this->matrix_free_handler->get_matrix_free());

  this->system_matrix->add_global_to_local_mapping(
    this->residual_global_to_local_solution);
  this->system_matrix->add_src_solution_subset(this->residual_src);

  this->update_system_matrix->add_global_to_local_mapping(
    this->newton_update_global_to_local_solution);
  this->update_system_matrix->add_src_solution_subset(this->newton_update_src);

  // Apply constraints
  this->constraint_handler->get_constraint(this->field_index)
    .distribute(*(this->solution_handler->get_solution_vector(this->field_index,
                                                              dependencyType::NORMAL)));
}

template <int dim, int degree>
inline void
identitySolver<dim, degree>::reinit()
{}

template <int dim, int degree>
inline void
identitySolver<dim, degree>::solve(const double &step_length)
{
  auto *solution = this->solution_handler->get_solution_vector(this->field_index,
                                                               dependencyType::NORMAL);

  // Compute the residual
  this->system_matrix->compute_residual(*this->residual, *solution);
  conditionalOStreams::pout_summary()
    << "  field: " << this->field_index
    << " Initial residual: " << this->residual->l2_norm() << std::flush;

  // Determine the residual tolerance
  this->compute_solver_tolerance();

  // Update solver controls
  this->solver_control.set_tolerance(this->tolerance);
  dealii::SolverCG<VectorType> cg_solver(this->solver_control);

  try
    {
      *this->newton_update = 0.0;
      cg_solver.solve(*(this->update_system_matrix),
                      *this->newton_update,
                      *this->residual,
                      dealii::PreconditionIdentity());
    }
  catch (...)
    {
      conditionalOStreams::pout_base()
        << "Warning: linear solver did not converge as per set tolerances.\n";
    }
  this->constraint_handler->get_constraint(this->field_index)
    .set_zero(*this->newton_update);

  conditionalOStreams::pout_summary()
    << " Final residual: " << this->solver_control.last_value()
    << " Steps: " << this->solver_control.last_step() << "\n"
    << std::flush;

  // Update the solutions
  (*solution).add(step_length, *this->newton_update);
  this->solution_handler->update(fieldSolveType::NONEXPLICIT_LINEAR, this->field_index);

  // Apply constraints
  // This may be redundant with the constraints on the update step.
  this->constraint_handler->get_constraint(this->field_index).distribute(*solution);
}

PRISMS_PF_END_NAMESPACE