// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/explicit_base.h>

#include <prismspf/config.h>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class handles the explicit solves of all constant fields
 */
template <unsigned int dim, unsigned int degree>
class ExplicitConstantSolver : public ExplicitBase<dim, degree>
{
public:
  /**
   * \brief Constructor.
   */
  explicit ExplicitConstantSolver(const SolverContext<dim, degree> &_solver_context);

  /**
   * \brief Destructor.
   */
  ~ExplicitConstantSolver() override = default;

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

PRISMS_PF_END_NAMESPACE