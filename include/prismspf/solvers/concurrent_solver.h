// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/solvers/solver_base.h>

#include <prismspf/config.h>

#include <functional>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class ConcurrentSolver : public SolverBase<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  ConcurrentSolver(const SolverContext<dim, degree> &_solver_context,
                   const FieldSolveType             &_field_solve_type,
                   Types::Index                      _solve_priority = 0)
    : SolverBase<dim, degree, number>(_solver_context,
                                      _field_solve_type,
                                      _solve_priority) {};

  /**
   * @brief Destructor.
   */
  ~ConcurrentSolver() override = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  ConcurrentSolver(const ConcurrentSolver &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  ConcurrentSolver &
  operator=(const ConcurrentSolver &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  ConcurrentSolver(ConcurrentSolver &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  ConcurrentSolver &
  operator=(ConcurrentSolver &&solver_base) noexcept = delete;

  /**
   * @brief Initialize the solver.
   */
  void
  init() override
  {
    // Call the base class init
    this->SolverBase<dim, degree, number>::init();

    // If the FieldSolveType is constant or the solver is empty we can just return early.
    if (this->get_field_solve_type() == FieldSolveType::ExplicitConstant ||
        this->solver_is_empty())
      {
        return;
      }

    // Create the MatrixFreeOperator
    system_matrix =
      std::make_unique<typename SolverBase<dim, degree, number>::SystemMatrixType>(
        this->get_subset_attributes(),
        this->get_pde_operator());

    // Apply constraints
    for (const auto &[index, variable] : this->get_subset_attributes())
      {
        this->get_solution_handler()
          .apply_constraints(index, this->get_constraint_handler().get_constraint(index));
      }

    // Set up the user-implemented equations and create the residual vectors
    system_matrix->clear();
    system_matrix->initialize(this->get_matrix_free_handler().get_matrix_free());

    // Resize the global to local solution vector
    global_to_local_solution.resize(
      this->get_subset_attributes().begin()->second.get_dependency_set_rhs().size());
    for (auto &vector : global_to_local_solution)
      {
        vector.resize(this->get_subset_attributes()
                        .begin()
                        ->second.get_dependency_set_rhs()
                        .begin()
                        ->size(),
                      Numbers::invalid_index);
      }

    // Create the subset of solution vectors and add the mapping to MatrixFreeOperator
    Types::Index dependency_index = 0;
    for (const auto &inner_dependency_set :
         this->get_subset_attributes().begin()->second.get_dependency_set_rhs())
      {
        Types::Index dependency_type = 0;
        for (const auto &field_type : inner_dependency_set)
          {
            // Skip if an invalid field type is found or the global_to_local_solution
            // already has an entry for this dependency index and dependency type
            if (field_type == Numbers::invalid_field_type ||
                global_to_local_solution[dependency_index][dependency_type] !=
                  Numbers::invalid_index)
              {
                dependency_type++;
                continue;
              }

            solution_subset.push_back(this->get_solution_handler().get_solution_vector(
              dependency_index,
              static_cast<DependencyType>(dependency_type)));
            new_solution_subset.push_back(
              this->get_solution_handler().get_new_solution_vector(dependency_index));
            global_to_local_solution[dependency_index][dependency_type] =
              solution_subset.size() - 1;

            dependency_type++;
          }

        dependency_index++;
      }
    system_matrix->add_global_to_local_mapping(global_to_local_solution);
  };

  /**
   * @brief Reinitialize the solver.
   */
  void
  reinit() override
  {
    // Call the base class reinit
    this->SolverBase<dim, degree, number>::reinit();

    // If the FieldSolveType is constant or the solver is empty we can just return early.
    if (this->get_field_solve_type() == FieldSolveType::ExplicitConstant ||
        this->solver_is_empty())
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

    // If the FieldSolveType is constant or the solver is empty we can just return early.
    if (this->get_field_solve_type() == FieldSolveType::ExplicitConstant ||
        this->solver_is_empty())
      {
        return;
      }
  };

  /**
   * @brief Print information about the solver to summary.log.
   */
  void
  print()
  {
    // Print the base class information
    this->SolverBase<dim, degree, number>::print();
  }

  /**
   * @brief Solve the explicit equations
   *
   * This is a common function for solving explicit (RHS only) equations that are
   * independent of one another.
   *
   * Rather than duplicate this code a bunch of times for explicit, postprocess, amr,
   * nucleation, etc... fields we have it here. Importantly, some of these have different
   * functions, so we require this function.
   */
  void
  solve_explicit_equations(
    const std::function<
      void(std::vector<typename SolverBase<dim, degree, number>::VectorType *> &,
           const std::vector<typename SolverBase<dim, degree, number>::VectorType *> &)>
      &function)
  {
    // Compute the update with the provided function
    function(new_solution_subset, solution_subset);

    // Scale the update by the respective (Scalar/Vector) invm. Note that we do this with
    // the original solution set to avoid some messy mapping.
    for (auto [index, vector] : this->get_solution_handler().get_new_solution_vector())
      {
        if (this->get_subset_attributes().find(index) !=
            this->get_subset_attributes().end())
          {
            vector->scale(this->get_invm_handler().get_invm(index));
          }
      }

    // Update the solutions
    this->get_solution_handler().update(this->get_field_solve_type());

    // Apply constraints
    // TODO (landinjm): This applies the constraints even to the old fields, which is
    // incorrect.
    for (const auto &[index, variable] : this->get_subset_attributes())
      {
        this->get_solution_handler()
          .apply_constraints(index, this->get_constraint_handler().get_constraint(index));
      }

    // Update the ghosts
    Timer::start_section("Update ghosts");
    this->get_solution_handler().update_ghosts();
    Timer::end_section("Update ghosts");
  }

  /**
   * @brief Get the system matrix.
   */
  [[nodiscard]] std::unique_ptr<
    typename SolverBase<dim, degree, number>::SystemMatrixType> &
  get_system_matrix()
  {
    return system_matrix;
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
   * @brief Matrix-free operator.
   */
  std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>
    system_matrix;

  /**
   * @brief Mapping from global solution vectors to the local ones
   */
  std::vector<std::vector<Types::Index>> global_to_local_solution;

  /**
   * @brief Subset of solutions fields that are necessary for concurrent solves.
   */
  std::vector<typename SolverBase<dim, degree, number>::VectorType *> solution_subset;

  /**
   * @brief Subset of new solutions fields that are necessary for concurrent solves.
   */
  std::vector<typename SolverBase<dim, degree, number>::VectorType *> new_solution_subset;
};

PRISMS_PF_END_NAMESPACE
