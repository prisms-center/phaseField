// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/solvers/linear_solver_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Class that handles the assembly and solving of a field with the identity
 * preconditioner (no preconditioner)
 */
template <unsigned int dim, unsigned int degree, typename number>
class IdentitySolver : public LinearSolverBase<dim, degree, number>
{
public:
  using SystemMatrixType = MatrixFreeOperator<dim, degree, number>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<number>;

  /**
   * @brief Constructor.
   */
  IdentitySolver(const UserInputParameters<dim>               &_user_inputs,
                 const VariableAttributes                     &_variable_attributes,
                 const MatrixfreeHandler<dim, number>         &_matrix_free_handler,
                 const ConstraintHandler<dim, degree, number> &_constraint_handler,
                 SolutionHandler<dim, number>                 &_solution_handler,
                 std::shared_ptr<const PDEOperator<dim, degree, number>> _pde_operator);

  /**
   * @brief Destructor.
   */
  ~IdentitySolver() override = default;

  /**
   * @brief Initialize the system.
   */
  void
  init() override;

  /**
   * @brief Reinitialize the system.
   */
  void
  reinit() override;

  /**
   * @brief Solve the system Ax=b.
   */
  void
  solve(const number &step_length = 1.0) override;
};

PRISMS_PF_END_NAMESPACE