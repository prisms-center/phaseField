// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/timer.h>

#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/linear_solver_identity.h>
#include <prismspf/solvers/solver_base.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <Types::Index dim, Types::Index degree, typename number>
class SequentialSolver : public SolverBase<dim, degree, number>
{
public:
  /**
   * @brief Constructor.
   */
  SequentialSolver(const SolverContext<dim, degree> &_solver_context,
                   const FieldSolveType             &_field_solve_type,
                   Types::Index                      _solve_priority = 0)
    : SolverBase<dim, degree, number>(_solver_context,
                                      _field_solve_type,
                                      _solve_priority) {};

  /**
   * @brief Destructor.
   */
  ~SequentialSolver() override = default;

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
  print() override
  {
    // Print the base class information
    this->SolverBase<dim, degree, number>::print();
  };

  /**
   * @brief Init a linear solver object of a given VariableAttributes.
   *
   * @param[in] variable The VariableAttributes
   */
  void
  init_linear_solver(const VariableAttributes &variable)
  {
    // Grab the global field index
    const Types::Index global_field_index = variable.get_field_index();

    if (this->get_user_inputs()
          .get_linear_solve_parameters()
          .get_linear_solve_parameters(global_field_index)
          .preconditioner == PreconditionerType::GMG)
      {
        gmg_solvers.emplace(
          global_field_index,
          std::make_unique<GMGSolver<dim, degree>>(this->get_user_inputs(),
                                                   variable,
                                                   this->get_matrix_free_handler(),
                                                   this->get_constraint_handler(),
                                                   this->get_triangulation_handler(),
                                                   this->get_dof_handler(),
                                                   this->get_mg_matrix_free_handler(),
                                                   this->get_solution_handler(),
                                                   this->get_pde_operator(),
                                                   this->get_pde_operator_float(),
                                                   this->get_mg_info()));
        gmg_solvers.at(global_field_index)->init();
      }
    else
      {
        identity_solvers.emplace(
          global_field_index,
          std::make_unique<IdentitySolver<dim, degree>>(this->get_user_inputs(),
                                                        variable,
                                                        this->get_matrix_free_handler(),
                                                        this->get_constraint_handler(),
                                                        this->get_solution_handler(),
                                                        this->get_pde_operator()));
        identity_solvers.at(global_field_index)->init();
      }
  }

  /**
   * @brief Init a explicit solver objects of a given VariableAttributes.
   *
   * @param[in] variable The VariableAttributes
   */
  void
  init_explicit_solver(const VariableAttributes &variable)
  {
    // Grab the global field index
    const Types::Index global_field_index = variable.get_field_index();

    // Creating temporary map to match types
    std::map<Types::Index, VariableAttributes> temp;
    temp.emplace(global_field_index, variable);
    subset_attributes_list.push_back(temp);

    // Create the implementation of MatrixFreeOperator with the subset of variable
    // attributes
    system_matrix[global_field_index] =
      std::make_unique<typename SolverBase<dim, degree, number>::SystemMatrixType>(
        subset_attributes_list.back(),
        this->get_pde_operator(),
        global_field_index);

    // Set up the user-implemented equations and create the residual vectors
    system_matrix[global_field_index]->clear();
    system_matrix[global_field_index]->initialize(
      this->get_matrix_free_handler().get_matrix_free());

    // Grab some data from the VariableAttributes
    const Types::Index max_fields =
      this->get_subset_attributes().begin()->second.get_max_fields();
    const Types::Index max_dependency_types =
      this->get_subset_attributes().begin()->second.get_max_dependency_types();

    // Resize the global to local solution vector
    global_to_local_solution[global_field_index].resize(max_fields * max_dependency_types,
                                                        Numbers::invalid_index);

    // Create the subset of solution vectors and add the mapping to MatrixFreeOperator
    new_solution_subset[global_field_index].push_back(
      this->get_solution_handler().get_new_solution_vector(global_field_index));
    solution_subset[global_field_index].push_back(
      this->get_solution_handler().get_solution_vector(global_field_index,
                                                       DependencyType::Normal));
    global_to_local_solution[global_field_index]
                            [global_field_index * max_dependency_types +
                             static_cast<Types::Index>(DependencyType::Normal)] = 0;

    Types::Index variable_index = 0;
    for (const auto &inner_dependency_set : variable.get_dependency_set_rhs())
      {
        Types::Index dependency_type = 0;
        for (const auto &field_type : inner_dependency_set)
          {
            // Skip if an invalid field type is found or the global_to_local_solution
            // already has an entry for this dependency index and dependency type
            if (field_type == Numbers::invalid_field_type ||
                global_to_local_solution[global_field_index]
                                        [variable_index * max_dependency_types +
                                         dependency_type] != Numbers::invalid_index)
              {
                dependency_type++;
                continue;
              }

            solution_subset[global_field_index].push_back(
              this->get_solution_handler().get_solution_vector(
                variable_index,
                static_cast<DependencyType>(dependency_type)));
            global_to_local_solution[global_field_index]
                                    [variable_index * max_dependency_types +
                                     dependency_type] =
                                      solution_subset.at(global_field_index).size() - 1;

            dependency_type++;
          }

        variable_index++;
      }
    system_matrix[global_field_index]->add_global_to_local_mapping(
      global_to_local_solution[global_field_index]);
  }

  /**
   * @brief Reinit a linear solver object of a given VariableAttributes.
   *
   * @param[in] variable The VariableAttributes
   */
  void
  reinit_linear_solver(const VariableAttributes &variable)
  {
    // Grab the global field index
    const Types::Index global_field_index = variable.get_field_index();

    if (this->get_user_inputs()
          .get_linear_solve_parameters()
          .get_linear_solve_parameters(global_field_index)
          .preconditioner == PreconditionerType::GMG)
      {
        gmg_solvers.at(global_field_index)->reinit();
      }
    else
      {
        identity_solvers.at(global_field_index)->reinit();
      }
  }

  /**
   * @brief Reinit a explicit solver objects of a given VariableAttributes.
   *
   * @param[in] variable The VariableAttributes
   */
  void
  reinit_explicit_solver([[maybe_unused]] const VariableAttributes &variable)
  {
    // Do nothing
  }

  /**
   * @brief Solve the explicit solver objects of a given VariableAttributes.
   *
   * @param[in] variable The VariableAttributes
   */
  void
  solve_explicit_solver(const VariableAttributes &variable)
  {
    // Grab the global field index
    const Types::Index global_field_index = variable.get_field_index();

    // Compute the update
    system_matrix[global_field_index]->compute_nonexplicit_auxiliary_update(
      new_solution_subset.at(global_field_index),
      solution_subset.at(global_field_index));

    // Scale the update by the respective (Scalar/Vector) invm.
    new_solution_subset.at(global_field_index)
      .at(0)
      ->scale(this->get_invm_handler().get_invm(global_field_index));

    // Update the solutions
    this->get_solution_handler().update(this->get_field_solve_type(), global_field_index);

    // Apply constraints
    this->get_constraint_handler()
      .get_constraint(global_field_index)
      .distribute(
        *(this->get_solution_handler().get_solution_vector(global_field_index,
                                                           DependencyType::Normal)));

    // Update the ghosts
    Timer::start_section("Update ghosts");
    this->get_solution_handler().update_ghosts();
    Timer::end_section("Update ghosts");
  }

  /**
   * @brief Solve the linear solver objects of a given VariableAttributes.
   *
   * @param[in] variable The VariableAttributes
   */
  void
  solve_linear_solver(const VariableAttributes &variable)
  {
    // Grab the global field index
    const Types::Index global_field_index = variable.get_field_index();

    // Skip if the field type is ImplicitTimeDependent and the increment is 0.
    if (variable.get_pde_type() == PDEType::ImplicitTimeDependent &&
        this->get_user_inputs().get_temporal_discretization().get_increment() == 0)
      {
        return;
      }

    if (this->get_user_inputs()
          .get_linear_solve_parameters()
          .get_linear_solve_parameters(global_field_index)
          .preconditioner == PreconditionerType::GMG)
      {
        gmg_solvers.at(global_field_index)->solve();
      }
    else
      {
        identity_solvers.at(global_field_index)->solve();
      }

    // Update the solutions
    this->get_solution_handler().update(this->get_field_solve_type(), global_field_index);

    // Update the ghosts
    Timer::start_section("Update ghosts");
    this->get_solution_handler().update_ghosts();
    Timer::end_section("Update ghosts");
  }

  /**
   * @brief Solve the linear solver objects of a given VariableAttributes.
   *
   * This function is overload specialized for co-nonlinear solves to use a given step
   * length and return a norm of the newton update. Additionally, it doesn't update the
   * solution.
   *
   * @param[in] variable The VariableAttributes
   * @param[in] step_length The step length of the linear solve. This is only used for
   * nonlinear solves when we don't want to use the entire solution.
   */
  double
  solve_linear_solver(const VariableAttributes &variable, const double &step_length)
  {
    // Grab the global field index
    const Types::Index global_field_index = variable.get_field_index();

    // Skip if the field type is ImplicitTimeDependent and the increment is 0.
    if (variable.get_pde_type() == PDEType::ImplicitTimeDependent &&
        this->get_user_inputs().get_temporal_discretization().get_increment() == 0)
      {
        return 0.0;
      }

    if (this->get_user_inputs()
          .get_linear_solve_parameters()
          .get_linear_solve_parameters(global_field_index)
          .preconditioner == PreconditionerType::GMG)
      {
        gmg_solvers.at(global_field_index)->solve(step_length);
      }
    else
      {
        identity_solvers.at(global_field_index)->solve(step_length);
      }

    // Update the ghosts
    Timer::start_section("Update ghosts");
    this->get_solution_handler().update_ghosts();
    Timer::end_section("Update ghosts");

    // Return the norm of the newton update
    if (this->get_user_inputs()
          .get_linear_solve_parameters()
          .get_linear_solve_parameters(global_field_index)
          .preconditioner == PreconditionerType::GMG)
      {
        return gmg_solvers.at(global_field_index)->get_newton_update_l2_norm();
      }
    return identity_solvers.at(global_field_index)->get_newton_update_l2_norm();
  }

  /**
   * @brief Get the matrix-free operator for the residual side.
   */
  [[nodiscard]] std::map<
    Types::Index,
    std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>> &
  get_system_matrix()
  {
    return system_matrix;
  }

  /**
   * @brief Get the matrix-free operator for the newton update side.
   */
  [[nodiscard]] std::map<
    Types::Index,
    std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>> &
  get_update_system_matrix()
  {
    return update_system_matrix;
  }

  /**
   * @brief Get the mapping from global solution vectors to the local ones.
   */
  [[nodiscard]] const std::map<Types::Index, std::vector<Types::Index>> &
  get_global_to_local_solution_mapping()
  {
    return global_to_local_solution;
  }

  /**
   * @brief Get the src solution subset.
   */
  [[nodiscard]] const std::map<
    Types::Index,
    std::vector<typename SolverBase<dim, degree, number>::VectorType *>> &
  get_src_solution_subset()
  {
    return solution_subset;
  }

  /**
   * @brief Get the dst solution subset.
   */
  [[nodiscard]] std::map<
    Types::Index,
    std::vector<typename SolverBase<dim, degree, number>::VectorType *>> &
  get_dst_solution_subset()
  {
    return new_solution_subset;
  }

private:
  /**
   * @brief Matrix-free operator for the residual side.
   */
  std::map<Types::Index,
           std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>>
    system_matrix;

  /**
   * @brief Matrix-free operator for the newton update side.
   */
  std::map<Types::Index,
           std::unique_ptr<typename SolverBase<dim, degree, number>::SystemMatrixType>>
    update_system_matrix;

  /**
   * @brief Mapping from global solution vectors to the local ones
   */
  std::map<Types::Index, std::vector<Types::Index>> global_to_local_solution;

  /**
   * @brief Subset of solutions fields that are necessary for concurrent solves.
   */
  std::map<Types::Index,
           std::vector<typename SolverBase<dim, degree, number>::VectorType *>>
    solution_subset;

  /**
   * @brief Subset of new solutions fields that are necessary for concurrent solves.
   */
  std::map<Types::Index,
           std::vector<typename SolverBase<dim, degree, number>::VectorType *>>
    new_solution_subset;

  /**
   * @brief List of subset attributes.
   */
  std::vector<std::map<Types::Index, VariableAttributes>> subset_attributes_list;

  /**
   * @brief Map of identity linear solvers
   */
  std::map<Types::Index, std::unique_ptr<IdentitySolver<dim, degree>>> identity_solvers;

  /**
   * @brief Map of geometric multigrid linear solvers
   */
  std::map<Types::Index, std::unique_ptr<GMGSolver<dim, degree>>> gmg_solvers;
};

PRISMS_PF_END_NAMESPACE
