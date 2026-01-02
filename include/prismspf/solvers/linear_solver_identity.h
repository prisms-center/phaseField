// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/solvers/linear_solver_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolverContext;

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
  IdentitySolver(const SolverContext<dim, degree, number> &_solver_context,
                 const VariableAttributes                 &_variable_attributes);

  /**
   * @brief Destructor.
   */
  ~IdentitySolver() override = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  IdentitySolver(const IdentitySolver &solver) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  IdentitySolver &
  operator=(const IdentitySolver &solver) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  IdentitySolver(IdentitySolver &&solver) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  IdentitySolver &
  operator=(IdentitySolver &&solver) noexcept = delete;

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
