// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/solver_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class SequentialSolver : public SolverBase<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  SequentialSolver(const SolverContext<dim, degree> &_solver_context,
                   const FieldSolveType             &_field_solve_type,
                   unsigned int                      _solve_priority = 0)
    : SolverBase<dim, degree, number>(_solver_context,
                                      _field_solve_type,
                                      _solve_priority) {};

  /**
   * @brief Destructor.
   */
  virtual ~SequentialSolver() = 0;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialSolver(const SequentialSolver &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  SequentialSolver &
  operator=(const SequentialSolver &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialSolver(SequentialSolver &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  SequentialSolver &
  operator=(SequentialSolver &&solver_base) noexcept = delete;

  /**
   * @brief Initialize the solver.
   */
  void
  init() override
  {
    // Call the base class init
    this->SolverBase<dim, degree, number>::init();

    // If the solver is empty we can just return early.
    if (this->solver_is_empty())
      {
        return;
      }
  };

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override
  {
    // Call the base class reinit
    this->SolverBase<dim, degree, number>::reinit();

    // If the solver is empty we can just return early.
    if (this->solver_is_empty())
      {
        return;
      }
  };

  /**
   * @brief Solve for a single update step.
   */
  void
  solve() override
  {
    // Call the base class solve
    this->SolverBase<dim, degree, number>::solve();

    // If the solver is empty we can just return early.
    if (this->solver_is_empty())
      {
        return;
      }
  };

  /**
   * @brief Print information about the solver to summary.log.
   */
  void
  print()
  { // Print the base class information
    this->SolverBase<dim, degree, number>::print();
  };

  /**
   * @brief Get the matrix-free operator for the residual side.
   */
  [[nodiscard]] std::map<
    unsigned int,
    std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>> &
  get_system_matrix()
  {
    return system_matrix;
  }

  /**
   * @brief Get the matrix-free operator for the newton update side.
   */
  [[nodiscard]] std::map<
    unsigned int,
    std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>> &
  get_update_system_matrix()
  {
    return update_system_matrix;
  }

  /**
   * @brief Get the mapping from global solution vectors to the local ones.
   */
  [[nodiscard]] const std::vector<std::vector<Types::Index>> &
  get_global_to_local_solution_mapping()
  {
    return global_to_local_solution;
  }

  /**
   * @brief Get the src solution subset.
   */
  [[nodiscard]] const std::vector<
    typename SolverBase<dim, degree, number>::VectorType *> &
  get_src_solution_subset()
  {
    return solution_subset;
  }

  /**
   * @brief Get the dst solution subset.
   */
  [[nodiscard]] std::vector<typename SolverBase<dim, degree, number>::VectorType *> &
  get_dst_solution_subset()
  {
    return new_solution_subset;
  }

private:
  /**
   * @brief Matrix-free operator for the residual side.
   */
  std::map<unsigned int,
           std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>>
    system_matrix;

  /**
   * @brief Matrix-free operator for the newton update side.
   */
  std::map<unsigned int,
           std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>>
    update_system_matrix;

  /**
   * @brief Mapping from global solution vectors to the local ones
   */
  std::map<unsigned int, std::vector<std::vector<Types::Index>>> global_to_local_solution;

  /**
   * @brief Subset of solutions fields that are necessary for concurrent solves.
   */
  std::map<unsigned int,
           std::vector<typename SolverBase<dim, degree, number>::VectorType *>>
    solution_subset;

  /**
   * @brief Subset of new solutions fields that are necessary for concurrent solves.
   */
  std::map<unsigned int,
           std::vector<typename SolverBase<dim, degree, number>::VectorType *>>
    new_solution_subset;

  /**
   * @brief List of subset attributes.
   */
  std::vector<std::map<unsigned int, VariableAttributes>> subset_attributes_list;

  /**
   * @brief Map of identity linear solvers
   */
  std::map<unsigned int, std::unique_ptr<IdentitySolver<dim, degree>>> identity_solvers;

  /**
   * @brief Map of geometric multigrid linear solvers
   */
  std::map<unsigned int, std::unique_ptr<GMGSolver<dim, degree>>> gmg_solvers;
};

PRISMS_PF_END_NAMESPACE
