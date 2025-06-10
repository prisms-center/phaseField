// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/nonexplicit_base.h>

#include <prismspf/config.h>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class handles all linear solves.
 */
template <unsigned int dim, unsigned int degree>
class NonexplicitLinearSolver : public NonexplicitBase<dim, degree>
{
public:
  using SystemMatrixType = MatrixFreeOperator<dim, degree, double>;

  /**
   * \brief Constructor.
   */
  NonexplicitLinearSolver(const SolverContext<dim, degree> &_solver_context,
                          const MGInfo<dim>                &_mg_info);

  /**
   * \brief Destructor.
   */
  ~NonexplicitLinearSolver() override = default;

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
  std::map<unsigned int, std::unique_ptr<IdentitySolver<dim, degree>>> identity_solvers;

  /**
   * \brief Map of geometric multigrid linear solvers
   */
  std::map<unsigned int, std::unique_ptr<GMGSolver<dim, degree>>> gmg_solvers;

  /**
   * \brief PDE operator but for floats!
   */
  std::shared_ptr<const PDEOperator<dim, degree, float>> pde_operator_float;

  /**
   * \brief Multigrid information
   */
  const MGInfo<dim> *mg_info;
};

PRISMS_PF_END_NAMESPACE
