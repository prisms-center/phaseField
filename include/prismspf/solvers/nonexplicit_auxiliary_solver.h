// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/nonexplicit_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief This class handles all auxiliary solves.
 */
template <unsigned int dim, unsigned int degree>
class NonexplicitAuxiliarySolver : public NonexplicitBase<dim, degree>
{
public:
  using SystemMatrixType = MatrixFreeOperator<dim, degree, double>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<double>;

  /**
   * \brief Constructor.
   */
  explicit NonexplicitAuxiliarySolver(const SolverContext<dim, degree> &_solver_context);

  /**
   * \brief Destructor.
   */
  ~NonexplicitAuxiliarySolver() override = default;

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
  std::map<unsigned int, std::map<std::pair<unsigned int, DependencyType>, unsigned int>>
    global_to_local_solution;

  /**
   * \brief Subset of solutions fields that are necessary for explicit solves.
   */
  std::map<unsigned int, std::vector<VectorType *>> solution_subset;

  /**
   * \brief Subset of new solutions fields that are necessary for explicit solves.
   */
  std::map<unsigned int, std::vector<VectorType *>> new_solution_subset;

  /**
   * \brief List of subset attributes.
   */
  std::vector<std::map<unsigned int, VariableAttributes>> subset_attributes_list;
};

PRISMS_PF_END_NAMESPACE
