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
template <int dim, int degree>
class explicitConstantSolver : public explicitBase<dim, degree>
{
public:
  /**
   * \brief Constructor.
   */
  explicitConstantSolver(
    const userInputParameters<dim>                         &_user_inputs,
    const matrixfreeHandler<dim, double>                   &_matrix_free_handler,
    const invmHandler<dim, degree, double>                 &_invm_handler,
    const constraintHandler<dim>                           &_constraint_handler,
    const dofHandler<dim>                                  &_dof_handler,
    const dealii::MappingQ1<dim>                           &_mapping,
    solutionHandler<dim>                                   &_solution_handler,
    std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator);

  /**
   * \brief Destructor.
   */
  ~explicitConstantSolver() override = default;

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