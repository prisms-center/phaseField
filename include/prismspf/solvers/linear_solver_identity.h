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

PRISMS_PF_END_NAMESPACE