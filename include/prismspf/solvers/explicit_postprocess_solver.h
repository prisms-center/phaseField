// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>

#include <prismspf/solvers/explicit_constant_solver.h>

#include <prismspf/config.h>

#ifdef PRISMS_PF_WITH_CALIPER
#  include <caliper/cali.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class handles the explicit solves of all postprocessed fields
 */
template <int dim, int degree>
class explicitPostprocessSolver : public explicitBase<dim, degree>
{
public:
  using SystemMatrixType = matrixFreeOperator<dim, degree, double>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<double>;

  /**
   * \brief Constructor.
   */
  explicitPostprocessSolver(
    const userInputParameters<dim>                         &_user_inputs,
    const matrixfreeHandler<dim>                           &_matrix_free_handler,
    const invmHandler<dim, degree>                         &_invm_handler,
    const constraintHandler<dim>                           &_constraint_handler,
    const dofHandler<dim>                                  &_dof_handler,
    const dealii::MappingQ1<dim>                           &_mapping,
    solutionHandler<dim>                                   &_solution_handler,
    std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator);

  /**
   * \brief Destructor.
   */
  ~explicitPostprocessSolver() = default;

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
   * \brief Mapping from global solution vectors to the local ones
   */
  std::unordered_map<std::pair<unsigned int, dependencyType>, unsigned int, pairHash>
    global_to_local_solution;

  /**
   * \brief Subset of solutions fields that are necessary for explicit solves.
   */
  std::vector<VectorType *> solution_subset;

  /**
   * \brief Subset of new solutions fields that are necessary for explicit solves.
   */
  std::vector<VectorType *> new_solution_subset;
};

PRISMS_PF_END_NAMESPACE