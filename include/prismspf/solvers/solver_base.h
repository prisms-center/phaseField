// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>

#include <map>
#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SolverBase
{
public:
  using SystemMatrixType = MatrixFreeOperator<dim, degree, number>;

  /**
   * \brief Constructor.
   */
  explicit SolverBase(const SolverContext<dim, degree> &_solver_context);

  /**
   * \brief Destructor.
   */
  virtual ~SolverBase() = 0;

  /**
   * \brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  SolverBase(const SolverBase &solver_base) = delete;

  /**
   * \brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  SolverBase &
  operator=(const SolverBase &solver_base) = delete;

  /**
   * \brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  SolverBase(SolverBase &&solver_base) noexcept = delete;

  /**
   * \brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  SolverBase &
  operator=(SolverBase &&solver_base) noexcept = delete;

  /**
   * \brief Initialize the solver.
   */
  virtual void
  init() = 0;

  /**
   * \brief Reinitialize the solver.
   */
  virtual void
  reinit() = 0;

  /**
   * \brief Solve for a single update step.
   */
  virtual void
  solve() = 0;

protected:
  /**
   * \brief Compute the subset of VariableAttributes that belongs to a given
   * FieldSolveType and solver order.
   *
   * This function creates a map of the VariablesAttributes that belong to a
   * FieldSolveType and solve order. The map can be accessed with get_subset_attributes.
   */
  void
  compute_subset_attributes(const FieldSolveType &field_solve_type,
                            unsigned int          solve_priority = 0);

private:
  /**
   * \brief Solver context.
   */
  const std::shared_ptr<SolverContext<dim, degree>> solver_context;

  /**
   * \brief Subset of variable attributes.
   */
  std::map<unsigned int, VariableAttributes> subset_attributes;

  /**
   * \brief Matrix-free operator
   */
  std::unique_ptr<SystemMatrixType> system_matrix;
};

PRISMS_PF_END_NAMESPACE
