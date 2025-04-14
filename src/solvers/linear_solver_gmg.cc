// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/quadrature_lib.h>
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
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/linear_solver_base.h>
#include <prismspf/solvers/linear_solver_gmg.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, int degree>
GMGSolver<dim, degree>::GMGSolver(
  const userInputParameters<dim>                         &_user_inputs,
  const variableAttributes                               &_variable_attributes,
  const matrixfreeHandler<dim>                           &_matrix_free_handler,
  const constraintHandler<dim>                           &_constraint_handler,
  const triangulationHandler<dim>                        &_triangulation_handler,
  const dofHandler<dim>                                  &_dof_handler,
  dealii::MGLevelObject<matrixfreeHandler<dim, float>>   &_mg_matrix_free_handler,
  solutionHandler<dim>                                   &_solution_handler,
  std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator,
  std::shared_ptr<const PDEOperator<dim, degree, float>>  _pde_operator_float)
  : linearSolverBase<dim, degree>(_user_inputs,
                                  _variable_attributes,
                                  _matrix_free_handler,
                                  _constraint_handler,
                                  _solution_handler,
                                  std::move(_pde_operator))
  , triangulation_handler(&_triangulation_handler)
  , dof_handler(&_dof_handler)
  , mg_matrix_free_handler(&_mg_matrix_free_handler)
  , pde_operator_float(std::move(_pde_operator_float))
{}

template <int dim, int degree>
GMGSolver<dim, degree>::~GMGSolver()
{
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      for (const auto *vector : mg_newton_update_src[level])
        {
          delete vector;
        }
    }
}

template <int dim, int degree>
inline void
GMGSolver<dim, degree>::init()
{
  // We solve the system with a global coarsening approach. There are two options when
  // doing this: geometric coarsening and polynomial coarsening. We only support geometric
  // as of now.

  // Grab the min and max level
  min_level = triangulation_handler->get_mg_min_level();
  max_level = triangulation_handler->get_mg_max_level();

  // Init the multilevel operator objects
  mg_operators =
    std::make_unique<dealii::MGLevelObject<LevelMatrixType>>(min_level,
                                                             max_level,
                                                             this->subset_attributes,
                                                             pde_operator_float,
                                                             this->field_index);
  mg_transfer_operators.resize(min_level, max_level);

  // Setup operator on each level
  mg_newton_update_src.resize(min_level, max_level);
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      // TODO (landinjm): Fix so mapping is same as rest of the problem. Do the same for
      // the finite element I think.
      // TODO (landinjm): This should include dof handlers for all dependency fields. That
      // also means I need some sort of local indexing.
      (*mg_matrix_free_handler)[level].reinit(
        mapping,
        dof_handler->get_mg_dof_handler(this->field_index, level),
        this->constraint_handler->get_mg_constraint(this->field_index, level),
        dealii::QGaussLobatto<1>(degree + 1));

      (*mg_operators)[level].initialize(
        (*mg_matrix_free_handler)[level].get_matrix_free());

      (*mg_operators)[level].add_global_to_local_mapping(
        this->newton_update_global_to_local_solution);

      // Setup src solutions for each level
      mg_newton_update_src[level].resize(this->newton_update_src.size());
      for (const auto &[pair, local_index] : this->newton_update_global_to_local_solution)
        {
          mg_newton_update_src[level][local_index] = new MGVectorType();
          (*mg_operators)[level]
            .initialize_dof_vector(*mg_newton_update_src[level][local_index], pair.first);
        }

      // Check that the vector partitioning is right
      // TODO (landinjm): Check this for all of the vectors. May also be able to remove
      // this
      Assert((*mg_operators)[level]
               .get_matrix_free()
               ->get_vector_partitioner(this->field_index)
               ->is_compatible(*(mg_newton_update_src[level][0]->get_partitioner())),
             dealii::ExcMessage("Incompatabile vector partitioners"));

      (*mg_operators)[level].add_src_solution_subset(mg_newton_update_src[level]);
    }
  mg_matrix = std::make_shared<dealii::mg::Matrix<MGVectorType>>(*mg_operators);

  // Setup transfer operators
  for (unsigned int level = min_level; level < max_level; ++level)
    {
      mg_transfer_operators[level + 1].reinit(
        dof_handler->get_mg_dof_handler(this->field_index, level + 1),
        dof_handler->get_mg_dof_handler(this->field_index, level),
        this->constraint_handler->get_mg_constraint(this->field_index, level + 1),
        this->constraint_handler->get_mg_constraint(this->field_index, level));
    }
  mg_transfer = std::make_shared<dealii::MGTransferGlobalCoarsening<dim, MGVectorType>>(
    mg_transfer_operators);

  this->system_matrix->clear();
  this->system_matrix->initialize(this->matrix_free_handler->get_matrix_free());
  this->update_system_matrix->clear();
  this->update_system_matrix->initialize(this->matrix_free_handler->get_matrix_free());

  this->system_matrix->add_global_to_local_mapping(
    this->residual_global_to_local_solution);
  this->system_matrix->add_src_solution_subset(this->residual_src);

  this->update_system_matrix->add_global_to_local_mapping(
    this->newton_update_global_to_local_solution);
  this->update_system_matrix->add_src_solution_subset(this->newton_update_src);

  // Apply constraints
  this->constraint_handler->get_constraint(this->field_index)
    .distribute(*(this->solution_handler->get_solution_vector(this->field_index,
                                                              dependencyType::NORMAL)));

  // TODO (landinjm): Should I put this somewhere else?
#ifdef DEBUG
  conditionalOStreams::pout_summary()
    << "\nMultigrid Setup Information for index " << this->field_index << ":\n"
    << "  Min level: " << min_level << "\n"
    << "  Max level: " << max_level << "\n";
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      conditionalOStreams::pout_summary()
        << "  Level: " << level << "\n"
        << "    Cells: "
        << triangulation_handler->get_mg_triangulation(level).n_global_active_cells()
        << "\n"
        << "    DoFs: "
        << dof_handler->get_mg_dof_handler(this->field_index, level).n_dofs() << "\n"
        << "    Constrained DoFs: "
        << this->constraint_handler->get_mg_constraint(this->field_index, level)
             .n_constraints()
        << "\n";
    }
  conditionalOStreams::pout_summary()
    << "  MG vertical communication efficiency: "
    << dealii::MGTools::vertical_communication_efficiency(
         triangulation_handler->get_mg_triangulation())
    << "\n"
    << "  MG workload imbalance: "
    << dealii::MGTools::workload_imbalance(triangulation_handler->get_mg_triangulation())
    << "\n\n"
    << std::flush;
#endif
}

template <int dim, int degree>
inline void
GMGSolver<dim, degree>::reinit()
{}

template <int dim, int degree>
inline void
GMGSolver<dim, degree>::solve(const double &step_length)
{
  const auto *current_dof_handler = dof_handler->get_dof_handlers().at(this->field_index);
  auto       *solution = this->solution_handler->get_solution_vector(this->field_index,
                                                               dependencyType::NORMAL);

  // Compute the residual
  this->system_matrix->compute_residual(*this->residual, *solution);
  conditionalOStreams::pout_summary()
    << "  field: " << this->field_index
    << " Initial residual: " << this->residual->l2_norm() << std::flush;

  // Determine the residual tolerance
  this->compute_solver_tolerance();

  // Update solver controls
  this->solver_control.set_tolerance(this->tolerance);
  dealii::SolverCG<VectorType> cg_solver(this->solver_control);

  // Interpolate the newton update src vector to each multigrid level
  for (const auto &[pair, local_index] : this->newton_update_global_to_local_solution)
    {
      // Create a temporary collection of the the dst pointers
      dealii::MGLevelObject<MGVectorType> mg_src_subset(min_level, max_level);
      for (unsigned int level = min_level; level < max_level; ++level)
        {
          mg_src_subset[level] = *mg_newton_update_src[level][local_index];
        }

      // Interpolate
      mg_transfer->interpolate_to_mg(*dof_handler->get_dof_handlers().at(pair.first),
                                     mg_src_subset,
                                     *this->newton_update_src[local_index]);
    }

  // Create smoother for each level
  using SmootherType = dealii::PreconditionChebyshev<LevelMatrixType, MGVectorType>;
  dealii::MGSmootherPrecondition<LevelMatrixType, SmootherType, MGVectorType> mg_smoother;
  dealii::MGLevelObject<typename SmootherType::AdditionalData> smoother_data(min_level,
                                                                             max_level);
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      smoother_data[level].smoothing_range =
        this->user_inputs->linear_solve_parameters.linear_solve.at(this->field_index)
          .smoothing_range;
      smoother_data[level].degree =
        this->user_inputs->linear_solve_parameters.linear_solve.at(this->field_index)
          .smoother_degree;
      smoother_data[level].eig_cg_n_iterations =
        this->user_inputs->linear_solve_parameters.linear_solve.at(this->field_index)
          .eig_cg_n_iterations;
      (*mg_operators)[level].compute_diagonal(this->field_index);
      smoother_data[level].preconditioner =
        (*mg_operators)[level].get_matrix_diagonal_inverse();
      smoother_data[level].constraints.copy_from(
        this->constraint_handler->get_mg_constraint(this->field_index, level));
    }
  mg_smoother.initialize(*mg_operators, smoother_data);

  dealii::MGCoarseGridApplySmoother<MGVectorType> mg_coarse;
  mg_coarse.initialize(mg_smoother);

  // Create multigrid object
  dealii::Multigrid<MGVectorType> multigrid(
    *mg_matrix,
    mg_coarse,
    *mg_transfer,
    mg_smoother,
    mg_smoother,
    min_level,
    max_level,
    dealii::Multigrid<MGVectorType>::Cycle::v_cycle);

  // Create the preconditioner
  const dealii::PreconditionMG<dim,
                               MGVectorType,
                               dealii::MGTransferGlobalCoarsening<dim, MGVectorType>>
    preconditioner(*current_dof_handler, multigrid, *mg_transfer);

  try
    {
      *this->newton_update = 0.0;
      cg_solver.solve(*(this->update_system_matrix),
                      *this->newton_update,
                      *this->residual,
                      preconditioner);
    }
  catch (...)
    {
      conditionalOStreams::pout_base()
        << "Warning: linear solver did not converge as per set tolerances.\n";
    }
  this->constraint_handler->get_constraint(this->field_index)
    .set_zero(*this->newton_update);

  conditionalOStreams::pout_summary()
    << " Final residual: " << this->solver_control.last_value()
    << " Steps: " << this->solver_control.last_step() << "\n"
    << std::flush;

  // Update the solutions
  (*solution).add(step_length, *this->newton_update);
  this->solution_handler->update(fieldSolveType::NONEXPLICIT_LINEAR, this->field_index);

  // Apply constraints
  // This may be redundant with the constraints on the update step.
  this->constraint_handler->get_constraint(this->field_index).distribute(*solution);
}

INSTANTIATE_BI_TEMPLATE(GMGSolver)

PRISMS_PF_END_NAMESPACE