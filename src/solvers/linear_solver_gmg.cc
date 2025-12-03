// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/multigrid.h>

#include <prismspf/core/conditional_ostreams.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/multigrid_info.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/solution_output.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/linear_solver_base.h>
#include <prismspf/solvers/linear_solver_gmg.h>
#include <prismspf/solvers/solver_context.h>

#include <prismspf/utilities/element_volume.h>

#include <prismspf/config.h>

#include <functional>
#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
GMGSolver<dim, degree, number>::GMGSolver(
  const SolverContext<dim, degree, number> &_solver_context,
  const VariableAttributes                 &_variable_attributes)
  : LinearSolverBase<dim, degree, number>(_solver_context, _variable_attributes)
  , min_level(_solver_context.get_mg_info().get_mg_min_level())
  , max_level(_solver_context.get_mg_info().get_mg_max_level())
{}

template <unsigned int dim, unsigned int degree, typename number>
void
GMGSolver<dim, degree, number>::init()
{
  // Call the base class init
  this->LinearSolverBase<dim, degree, number>::init();

  // Basic intialization that is the same as the identity solve.
  this->clear_system_matrices();
  this->initialize_system_matrices();
  this->finalize_system_matrices();

  // Apply constraints
  this->apply_constraints();

  // We solve the system with a global coarsening approach. There are two options when
  // doing this: geometric coarsening and polynomial coarsening. We only support geometric
  // as of now.

  // Grab some data from the VariableAttributes
  const Types::Index max_fields = this->get_variable_attributes().get_max_fields();
  const Types::Index max_dependency_types =
    this->get_variable_attributes().get_max_dependency_types();

  // Print the local indices and global ones and check that they match the MGInfo provided
  // ones.
  for (Types::Index _field_index = 0; _field_index < max_fields; _field_index++)
    {
      for (Types::Index dependency_type = 0; dependency_type < max_dependency_types;
           dependency_type++)
        {
          Types::Index local_index =
            this->get_newton_update_global_to_local_solution()[(_field_index *
                                                                max_dependency_types) +
                                                               dependency_type];
          // Skip if the local index is invalid
          if (local_index == Numbers::invalid_index)
            {
              continue;
            }
          // Grab the local index for the change variable
          if (dependency_type == static_cast<Types::Index>(DependencyType::Change))
            {
              change_local_index = local_index;
            }
          for (unsigned int level = min_level; level <= max_level; ++level)
            {
              Assert(local_index ==
                       this->get_mg_info().get_local_index(_field_index, level),
                     dealii::ExcMessage(
                       "The multigrid info indexing must match the local to "
                       "global ones in this solver. "));
            }
        }
    }

  // Init the multilevel operator objects
  mg_operators = std::make_unique<dealii::MGLevelObject<LevelMatrixType>>(
    min_level,
    max_level,
    this->get_subset_attributes(),
    this->get_pde_operator_float(),
    this->get_variable_attributes().get_solve_block(),
    this->get_field_index(),
    true);

  // Setup operator on each level
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      (*mg_operators)[level].initialize(
        this->get_matrix_free_container().get_mg_matrix_free(level),
        this->get_element_volume_container().get_mg_element_volume(level),
        {change_local_index});

      (*mg_operators)[level].add_global_to_local_mapping(
        this->get_newton_update_global_to_local_solution());

      (*mg_operators)[level].add_src_solution_subset(
        this->get_solution_handler().get_mg_solution_vector(level));
    }
  mg_matrix = std::make_shared<dealii::mg::Matrix<MGVectorType>>(*mg_operators);

  // Setup transfer operators
  // For now I'll just make a bunch of them based on the local indices.
  // TODO (landinjm): This is awful please fix.
  mg_transfer_operators.resize(this->get_newton_update_global_to_local_solution().size());
  mg_transfer.resize(mg_transfer_operators.size());
  for (Types::Index _field_index = 0; _field_index < max_fields; _field_index++)
    {
      for (Types::Index dependency_type = 0; dependency_type < max_dependency_types;
           dependency_type++)
        {
          Types::Index local_index =
            this->get_newton_update_global_to_local_solution()[(_field_index *
                                                                max_dependency_types) +
                                                               dependency_type];
          // Skip if the local index is invalid
          if (local_index == Numbers::invalid_index)
            {
              continue;
            }

          mg_transfer_operators[local_index].resize(min_level, max_level);

          for (unsigned int level = min_level; level < max_level; ++level)
            {
              mg_transfer_operators[local_index][level + 1].reinit(
                (*this->get_dof_handler().get_mg_dof_handlers(level + 1)[local_index]),
                *(this->get_dof_handler().get_mg_dof_handlers(level)[local_index]),
                this->get_constraint_handler().get_mg_constraint(level + 1, local_index),
                this->get_constraint_handler().get_mg_constraint(level, local_index));
            }
          mg_transfer[local_index] =
            std::make_shared<dealii::MGTransferGlobalCoarsening<dim, MGVectorType>>(
              mg_transfer_operators[local_index],
              std::function<void(const unsigned int, MGVectorType &)>(
                [&](const unsigned int level, MGVectorType &vec)
                {
                  (*mg_operators)[level].initialize_dof_vector(vec, local_index);
                }));
        }
    }

  ConditionalOStreams::pout_summary()
    << "\nMultigrid Setup Information for index " << this->get_field_index() << ":\n"
    << "  Min level: " << min_level << "\n"
    << "  Max level: " << max_level << "\n"
    << "  MG vertical communication efficiency: "
    << dealii::MGTools::vertical_communication_efficiency(
         this->get_triangulation_handler().get_mg_triangulation())
    << "\n"
    << "  MG workload imbalance: "
    << dealii::MGTools::workload_imbalance(
         this->get_triangulation_handler().get_mg_triangulation())
    << "\n\n"
    << std::flush;
}

template <unsigned int dim, unsigned int degree, typename number>
void
GMGSolver<dim, degree, number>::reinit()
{
  // Call the base class reinit
  this->LinearSolverBase<dim, degree, number>::reinit();

  // Basic intialization that is the same as the identity solve.
  this->clear_system_matrices();
  this->initialize_system_matrices();
  this->finalize_system_matrices();

  // Apply constraints
  this->apply_constraints();

  // We solve the system with a global coarsening approach. There are two options when
  // doing this: geometric coarsening and polynomial coarsening. We only support geometric
  // as of now.

  // Grab some data from the VariableAttributes
  const Types::Index max_fields = this->get_variable_attributes().get_max_fields();
  const Types::Index max_dependency_types =
    this->get_variable_attributes().get_max_dependency_types();

  // Setup operator on each level
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      (*mg_operators)[level].clear();
      (*mg_operators)[level].initialize(
        this->get_matrix_free_container().get_mg_matrix_free(level),
        this->get_element_volume_container().get_mg_element_volume(level),
        {change_local_index});
      (*mg_operators)[level].add_global_to_local_mapping(
        this->get_newton_update_global_to_local_solution());
      (*mg_operators)[level].add_src_solution_subset(
        this->get_solution_handler().get_mg_solution_vector(level));
    }
  mg_matrix = std::make_shared<dealii::mg::Matrix<MGVectorType>>(*mg_operators);

  // Setup transfer operators
  // For now I'll just make a bunch of them based on the local indices.
  // TODO (landinjm): This is awful please fix.
  mg_transfer_operators.clear();
  mg_transfer.clear();
  mg_transfer_operators.resize(this->get_newton_update_global_to_local_solution().size());
  mg_transfer.resize(mg_transfer_operators.size());
  for (Types::Index _field_index = 0; _field_index < max_fields; _field_index++)
    {
      for (Types::Index dependency_type = 0; dependency_type < max_dependency_types;
           dependency_type++)
        {
          Types::Index local_index =
            this->get_newton_update_global_to_local_solution()[(_field_index *
                                                                max_dependency_types) +
                                                               dependency_type];
          // Skip if the local index is invalid
          if (local_index == Numbers::invalid_index)
            {
              continue;
            }

          mg_transfer_operators[local_index].resize(min_level, max_level);

          for (unsigned int level = min_level; level < max_level; ++level)
            {
              mg_transfer_operators[local_index][level + 1].reinit(
                *(this->get_dof_handler().get_mg_dof_handlers(level + 1)[local_index]),
                *(this->get_dof_handler().get_mg_dof_handlers(level)[local_index]),
                this->get_constraint_handler().get_mg_constraint(level + 1, local_index),
                this->get_constraint_handler().get_mg_constraint(level, local_index));
            }
          mg_transfer[local_index] =
            std::make_shared<dealii::MGTransferGlobalCoarsening<dim, MGVectorType>>(
              mg_transfer_operators[local_index],
              std::function<void(const unsigned int, MGVectorType &)>(
                [&](const unsigned int level, MGVectorType &vec)
                {
                  (*mg_operators)[level].initialize_dof_vector(vec, local_index);
                }));
        }
    }

  ConditionalOStreams::pout_summary()
    << "\nMultigrid Setup Information for index " << this->get_field_index() << ":\n"
    << "  Min level: " << min_level << "\n"
    << "  Max level: " << max_level << "\n"
    << "  MG vertical communication efficiency: "
    << dealii::MGTools::vertical_communication_efficiency(
         this->get_triangulation_handler().get_mg_triangulation())
    << "\n"
    << "  MG workload imbalance: "
    << dealii::MGTools::workload_imbalance(
         this->get_triangulation_handler().get_mg_triangulation())
    << "\n\n"
    << std::flush;
}

template <unsigned int dim, unsigned int degree, typename number>
void
GMGSolver<dim, degree, number>::solve(const number &step_length)
{
  const auto *current_dof_handler =
    this->get_dof_handler().get_dof_handlers().at(this->get_field_index());
  auto *solution =
    this->get_solution_handler().get_solution_vector(this->get_field_index(),
                                                     DependencyType::Normal);

  // Compute the residual
  this->get_system_matrix()->compute_residual(*this->get_residual(), *solution);
  if (this->get_user_inputs().get_output_parameters().should_output(
        this->get_user_inputs().get_temporal_discretization().get_increment()))
    {
      ConditionalOStreams::pout_summary()
        << "  field: " << this->get_field_index()
        << " Initial residual: " << this->get_residual()->l2_norm() << std::flush;
    }

  // Determine the residual tolerance
  this->compute_solver_tolerance();

  // Update solver controls
  this->get_solver_control().set_tolerance(static_cast<double>(this->get_tolerance()));
  dealii::SolverCG<VectorType> cg_solver(this->get_solver_control());

  // Grab some data from the VariableAttributes
  const Types::Index max_fields = this->get_variable_attributes().get_max_fields();
  const Types::Index max_dependency_types =
    this->get_variable_attributes().get_max_dependency_types();

  // Interpolate the newton update src vector to each multigrid level
  for (Types::Index _field_index = 0; _field_index < max_fields; _field_index++)
    {
      for (Types::Index dependency_type = 0; dependency_type < max_dependency_types;
           dependency_type++)
        {
          Types::Index local_index =
            this->get_newton_update_global_to_local_solution()[(_field_index *
                                                                max_dependency_types) +
                                                               dependency_type];
          // Skip if the local index is invalid
          if (local_index == Numbers::invalid_index)
            {
              continue;
            }

          // Create a temporary collection of the the dst pointers
          dealii::MGLevelObject<MGVectorType> mg_src_subset(min_level, max_level);
          for (unsigned int level = min_level; level <= max_level; ++level)
            {
              mg_src_subset[level] =
                *this->get_solution_handler().get_mg_solution_vector(level, local_index);

              // Check that the vector partitioning is right
              Assert((*mg_operators)[level]
                       .get_matrix_free()
                       ->get_vector_partitioner(local_index)
                       ->is_compatible(*mg_src_subset[level].get_partitioner()),
                     dealii::ExcMessage("Incompatabile vector partitioners"));
            }

          // Interpolate
          mg_transfer[local_index]->interpolate_to_mg(
            *(this->get_dof_handler().get_dof_handlers().at(_field_index)),
            mg_src_subset,
            *this->get_newton_update_src()[local_index]);

          // Copy back the vectors
          for (unsigned int level = min_level; level <= max_level; ++level)
            {
              *this->get_solution_handler().get_mg_solution_vector(level, local_index) =
                mg_src_subset[level];
            }
        }
    }

  // Create smoother for each level
  using SmootherType = dealii::PreconditionChebyshev<LevelMatrixType, MGVectorType>;
  dealii::MGSmootherPrecondition<LevelMatrixType, SmootherType, MGVectorType> mg_smoother;
  dealii::MGLevelObject<typename SmootherType::AdditionalData> smoother_data(min_level,
                                                                             max_level);
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      smoother_data[level].smoothing_range =
        this->get_user_inputs()
          .get_linear_solve_parameters()
          .get_linear_solve_parameters(this->get_field_index())
          .smoothing_range;
      smoother_data[level].degree =
        this->get_user_inputs()
          .get_linear_solve_parameters()
          .get_linear_solve_parameters(this->get_field_index())
          .smoother_degree;
      smoother_data[level].eig_cg_n_iterations =
        this->get_user_inputs()
          .get_linear_solve_parameters()
          .get_linear_solve_parameters(this->get_field_index())
          .eig_cg_n_iterations;
      (*mg_operators)[level].compute_diagonal(change_local_index);
      smoother_data[level].preconditioner =
        (*mg_operators)[level].get_matrix_diagonal_inverse();
      smoother_data[level].constraints.copy_from(
        this->get_constraint_handler().get_mg_constraint(level, change_local_index));
    }
  mg_smoother.initialize(*mg_operators, smoother_data);

  dealii::MGCoarseGridApplySmoother<MGVectorType> mg_coarse;
  mg_coarse.initialize(mg_smoother);

  // Create multigrid object
  dealii::Multigrid<MGVectorType> multigrid(
    *mg_matrix,
    mg_coarse,
    *mg_transfer[change_local_index],
    mg_smoother,
    mg_smoother,
    min_level,
    max_level,
    dealii::Multigrid<MGVectorType>::Cycle::v_cycle);

  // Create the preconditioner
  const dealii::PreconditionMG<dim,
                               MGVectorType,
                               dealii::MGTransferGlobalCoarsening<dim, MGVectorType>>
    preconditioner(*current_dof_handler, multigrid, *mg_transfer[change_local_index]);

  try
    {
      *this->get_newton_update() = 0.0;
      cg_solver.solve(*(this->get_update_system_matrix()),
                      *(this->get_newton_update()),
                      *(this->get_residual()),
                      preconditioner);
    }
  catch (...)
    {
      ConditionalOStreams::pout_base()
        << "Warning: linear solver did not converge as per set tolerances.\n";
    }
  this->get_constraint_handler()
    .get_constraint(this->get_field_index())
    .set_zero(*this->get_newton_update());

  if (this->get_user_inputs().get_output_parameters().should_output(
        this->get_user_inputs().get_temporal_discretization().get_increment()))
    {
      ConditionalOStreams::pout_summary()
        << " Final residual: " << this->get_solver_control().last_value()
        << " Steps: " << this->get_solver_control().last_step() << "\n"
        << std::flush;
    }

  // Update the solutions
  (*solution).add(step_length, *this->get_newton_update());

  // Apply constraints
  // This may be redundant with the constraints on the update step.
  this->get_constraint_handler()
    .get_constraint(this->get_field_index())
    .distribute(*solution);
}

#include "solvers/linear_solver_gmg.inst"

PRISMS_PF_END_NAMESPACE
