// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/explicit_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief This class handles the explicit solves of all explicit fields
 */
template <unsigned int dim, unsigned int degree>
class ExplicitSolver : public ExplicitBase<dim, degree>
{
public:
  using SystemMatrixType = MatrixFreeOperator<dim, degree, double>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<double>;

  /**
   * @brief Constructor.
   */
  explicit ExplicitSolver(const SolverContext<dim, degree> &_solver_context);

  /**
   * @brief Destructor.
   */
  ~ExplicitSolver() override = default;

  /**
   * @brief Initialize system.
   */
  void
  init() override;

  /**
   * @brief Solve a single update step.
   */
  void
  solve() override;

private:
  /**
   * @brief Mapping from global solution vectors to the local ones
   */
  std::vector<std::vector<Types::Index>> global_to_local_solution;

  /**
   * @brief Subset of solutions fields that are necessary for explicit solves.
   */
  std::vector<VectorType *> solution_subset;

  /**
   * @brief Subset of new solutions fields that are necessary for explicit solves.
   */
  std::vector<VectorType *> new_solution_subset;
};

PRISMS_PF_END_NAMESPACE
