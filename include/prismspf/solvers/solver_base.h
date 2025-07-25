// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/base/exceptions.h>
#include <deal.II/numerics/vector_tools.h>

#include <prismspf/core/initial_conditions.h>
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
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<number>;

  /**
   * @brief Constructor.
   */
  SolverBase(const SolverContext<dim, degree> &_solver_context,
             const FieldSolveType             &_field_solve_type,
             Types::Index                      _solve_priority = 0)
    : solver_context(std::make_shared<SolverContext<dim, degree>>(_solver_context))
    , field_solve_type(_field_solve_type)
    , solve_priority(_solve_priority)
  {}

  /**
   * @brief Destructor.
   */
  virtual ~SolverBase() = default;

  /**
   * @brief Copy constructor.
   *
   * Deleted so solver instances aren't copied.
   */
  SolverBase(const SolverBase &solver_base) = delete;

  /**
   * @brief Copy assignment.
   *
   * Deleted so solver instances aren't copied.
   */
  SolverBase &
  operator=(const SolverBase &solver_base) = delete;

  /**
   * @brief Move constructor.
   *
   * Deleted so solver instances aren't moved.
   */
  SolverBase(SolverBase &&solver_base) noexcept = delete;

  /**
   * @brief Move assignment.
   *
   * Deleted so solver instances aren't moved.
   */
  SolverBase &
  operator=(SolverBase &&solver_base) noexcept = delete;

  /**
   * @brief Initialize the solver.
   */
  virtual void
  init()
  {
    // Update the subset of variable attributes
    update_subset_attributes(field_solve_type, solve_priority);

    // If the subset attribute is empty return early
    if (solver_is_empty())
      {
        ConditionalOStreams::pout_base() << "  no fields for this solver exist\n"
                                         << std::flush;
        return;
      }

    // Set the initial condition
    set_initial_condition();

    // Apply constraints. This part is neccessary so they are taken into account for
    // adaptive meshing
    for (const auto &[index, variable] : subset_attributes)
      {
        get_constraint_handler().get_constraint(index).distribute(
          *(get_solution_handler().get_solution_vector(index, DependencyType::Normal)));
      }
  };

  /**
   * @brief Reinitialize the solver.
   */
  virtual void
  reinit()
  {
    // If the subset attribute is empty return early
    if (solver_is_empty())
      {
        return;
      }

    // Apply constraints. This part is neccessary so they are taken into account for
    // adaptive meshing
    for (const auto &[index, variable] : subset_attributes)
      {
        get_constraint_handler().get_constraint(index).distribute(
          *(get_solution_handler().get_solution_vector(index, DependencyType::Normal)));
      }
  };

  /**
   * @brief Solve for a single update step.
   */
  virtual void
  solve()
  {
    // If the subset attribute is empty return early
    if (solver_is_empty())
      {
        return;
      }

    // Do nothing
  };

  /**
   * @brief Print information about the solver to summary.log.
   */
  virtual void
  print();

  /**
   * @brief Whether the subset attributes is empty.
   *
   * This function is used so we can return early.
   */
  [[nodiscard]] bool
  solver_is_empty() const
  {
    return subset_attributes.empty();
  }

  /**
   * @brief Compute the subset of VariableAttributes that belongs to a given
   * FieldSolveType and solver order.
   *
   * This function creates and returns a map of the VariablesAttributes that belong to a
   * FieldSolveType and solve order.
   */
  [[nodiscard]] std::map<Types::Index, VariableAttributes>
  compute_subset_attributes(const FieldSolveType &field_solve_type,
                            Types::Index          solve_priority) const
  {
    std::map<Types::Index, VariableAttributes> local_subset_attributes;

    // TODO (landinjm): Use the solve priority
    (void) solve_priority;

    for (const auto &[index, variable] :
         solver_context->get_user_inputs().get_variable_attributes())
      {
        if (variable.get_field_solve_type() == field_solve_type)
          {
            local_subset_attributes.emplace(index, variable);
          }
      }

    return local_subset_attributes;
  };

  /**
   * @brief Compute and update the subset of VariableAttributes that belongs to a given
   * FieldSolveType and solver order.
   *
   * This function creates a map of the VariablesAttributes that belong to a
   * FieldSolveType and solve order. The map can be accessed with get_subset_attributes.
   */
  void
  update_subset_attributes(const FieldSolveType &field_solve_type,
                           Types::Index          solve_priority)
  {
    subset_attributes = compute_subset_attributes(field_solve_type, solve_priority);
  };

  /**
   * @brief Set the initial condition according to subset_attributes.
   *
   * This only sets the initial conditions for ExplicitTimeDependent,
   * ImplicitTimeDependent, TimeIndependent, and Constant fields.
   */
  void
  set_initial_condition()
  {
    for (const auto &[index, variable] : subset_attributes)
      {
        // TODO (landinjm): Skip certain fields for initial conditions

        Assert(solver_context->get_dof_handler().get_dof_handlers().size() > index,
               dealii::ExcMessage(
                 "The const DoFHandler set is smaller than the given index = " +
                 std::to_string(index)));
        Assert(subset_attributes.contains(index),
               dealii::ExcMessage(
                 "There is no entry in the attribute subset for the given index = " +
                 std::to_string(index)));

        if (solver_context->get_user_inputs()
              .get_load_initial_condition_parameters()
              .get_read_initial_conditions_from_file())
          {
            auto &initial_condition_parameters =
              solver_context->get_user_inputs().get_load_initial_condition_parameters();
            for (const auto &initial_condition_file :
                 initial_condition_parameters.get_initial_condition_files())
              {
                auto iterator =
                  std::find(initial_condition_file.simulation_variable_names.begin(),
                            initial_condition_file.simulation_variable_names.end(),
                            variable.get_name());
                if (iterator != initial_condition_file.simulation_variable_names.end())
                  {
                    dealii::VectorTools::interpolate(
                      solver_context->get_mapping(),
                      *(solver_context->get_dof_handler().get_dof_handlers().at(index)),
                      ReadInitialCondition<dim>(
                        initial_condition_file.filename + "." +
                          initial_condition_file.file_extension,
                        initial_condition_file.file_variable_names
                          [iterator -
                           initial_condition_file.simulation_variable_names.begin()],
                        subset_attributes.at(index).get_field_type()),
                      *(solver_context->get_solution_handler()
                          .get_solution_vector(index, DependencyType::Normal)));
                  }
              }
          }
        else
          {
            dealii::VectorTools::interpolate(
              solver_context->get_mapping(),
              *(solver_context->get_dof_handler().get_dof_handlers().at(index)),
              InitialCondition<dim, degree>(index,
                                            subset_attributes.at(index).get_field_type(),
                                            solver_context->get_pde_operator()),
              *(solver_context->get_solution_handler()
                  .get_solution_vector(index, DependencyType::Normal)));
          }

        // TODO (landinjm): Fix so that we apply some sort of initial condition to all old
        // vector for all types.
        solver_context->get_solution_handler().apply_initial_condition_for_old_fields();
      }
  };

  /**
   * @brief Get the user-inputs.
   */
  [[nodiscard]] const UserInputParameters<dim> &
  get_user_inputs() const
  {
    return solver_context->get_user_inputs();
  }

  /**
   * @brief Get the matrix-free object handler for non-multigrid data.
   */
  [[nodiscard]] const MatrixfreeHandler<dim, double> &
  get_matrix_free_handler() const
  {
    return solver_context->get_matrix_free_handler();
  }

  /**
   * @brief Get the triangulation handler.
   */
  [[nodiscard]] const TriangulationHandler<dim> &
  get_triangulation_handler() const
  {
    return solver_context->get_triangulation_handler();
  }

  /**
   * @brief Get the invm handler.
   */
  [[nodiscard]] const InvmHandler<dim, degree, double> &
  get_invm_handler() const
  {
    return solver_context->get_invm_handler();
  }

  /**
   * @brief Get the constraint handler.
   */
  [[nodiscard]] const ConstraintHandler<dim, degree> &
  get_constraint_handler() const
  {
    return solver_context->get_constraint_handler();
  }

  /**
   * @brief Get the dof handler.
   */
  [[nodiscard]] const DofHandler<dim> &
  get_dof_handler() const
  {
    return solver_context->get_dof_handler();
  }

  /**
   * @brief Get the mapping.
   */
  [[nodiscard]] const dealii::MappingQ1<dim> &
  get_mapping() const
  {
    return solver_context->get_mapping();
  }

  /**
   * @brief Get the multigrid info.
   */
  [[nodiscard]] const MGInfo<dim> &
  get_mg_info() const
  {
    return solver_context->get_mg_info();
  }

  /**
   * @brief Get the mg matrix-free handler.
   */
  [[nodiscard]] dealii::MGLevelObject<MatrixfreeHandler<dim, float>> &
  get_mg_matrix_free_handler()
  {
    return solver_context->get_mg_matrix_free_handler();
  }

  /**
   * @brief Get the solution handler.
   */
  [[nodiscard]] SolutionHandler<dim> &
  get_solution_handler() const
  {
    return solver_context->get_solution_handler();
  }

  /**
   * @brief Get the subset attributes.
   */
  [[nodiscard]] const std::map<Types::Index, VariableAttributes> &
  get_subset_attributes() const
  {
    return subset_attributes;
  }

  /**
   * @brief Get the pde operator.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperator<dim, degree, double>> &
  get_pde_operator() const
  {
    return solver_context->get_pde_operator();
  }

  /**
   * @brief Get the pde operator for float.
   */
  [[nodiscard]] const std::shared_ptr<const PDEOperator<dim, degree, float>> &
  get_pde_operator_float() const
  {
    return solver_context->get_pde_operator_float();
  }

  /**
   * @brief Get the field solve type.
   */
  [[nodiscard]] const FieldSolveType &
  get_field_solve_type() const
  {
    Assert(field_solve_type != Numbers::invalid_field_solve_type,
           dealii::ExcMessage("The field solve type is invalid."));
    return field_solve_type;
  }

private:
  /**
   * @brief Solver context.
   */
  const std::shared_ptr<SolverContext<dim, degree>> solver_context;

  /**
   * @brief Field solve type.
   */
  const FieldSolveType field_solve_type = Numbers::invalid_field_solve_type;

  /**
   * @brief Solve priority.
   */
  const Types::Index solve_priority = Numbers::invalid_index;

  /**
   * @brief Subset of variable attributes.
   */
  std::map<Types::Index, VariableAttributes> subset_attributes;
};

PRISMS_PF_END_NAMESPACE
