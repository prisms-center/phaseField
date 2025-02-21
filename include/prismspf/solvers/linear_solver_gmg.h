// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef linear_solver_gmg_h
#define linear_solver_gmg_h

#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_global_coarsening.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <prismspf/config.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/solution_output.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/solvers/linear_solver_base.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * Forward declaration for user-implemented PDE class.
 */
template <int dim, int degree, typename number>
class customPDE;

/**
 * \brief Class that handles the assembly and solving of a field with a GMG preconditioner
 */
template <int dim, int degree>
class GMGSolver : public linearSolverBase<dim, degree>
{
public:
  using SystemMatrixType = customPDE<dim, degree, double>;
  using LevelMatrixType  = customPDE<dim, degree, float>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<double>;
  using MGVectorType     = dealii::LinearAlgebra::distributed::Vector<float>;

  /**
   * \brief Constructor.
   */
  GMGSolver(const userInputParameters<dim>                       &_user_inputs,
            const variableAttributes                             &_variable_attributes,
            const matrixfreeHandler<dim>                         &_matrix_free_handler,
            const constraintHandler<dim>                         &_constraint_handler,
            const triangulationHandler<dim>                      &_triangulation_handler,
            const prisms::dofHandler<dim>                        &_dof_handler,
            dealii::MGLevelObject<matrixfreeHandler<dim, float>> &_mg_matrix_free_handler,
            solutionHandler<dim>                                 &_solution_handler);

  /**
   * \brief Destructor.
   */
  ~GMGSolver() override;

  /**
   * \brief Initialize the system.
   */
  void
  init() override;

  /**
   * \brief Reinitialize the system.
   */
  void
  reinit() override;

  /**
   * \brief Solve the system Ax=b.
   */
  void
  solve(const double step_length = 1.0) override;

private:
  /**
   * \brief Triangulation handler.
   */
  const triangulationHandler<dim> &triangulation_handler;

  /**
   * \brief DoF handler.
   */
  const prisms::dofHandler<dim> &dof_handler;

  /**
   * \brief Matrix-free object handler for multigrid data.
   */
  dealii::MGLevelObject<matrixfreeHandler<dim, float>> &mg_matrix_free_handler;

  /**
   * \brief Minimum multigrid level
   */
  unsigned int min_level = 0;

  /**
   * \brief Maximum multigrid level
   */
  unsigned int max_level = 0;

  /**
   * \brief Mappings to and from reference cell.
   */
  const dealii::MappingQ1<dim> mapping;

  /**
   * \brief Collection of triangulations for each multigrid level.
   */
  std::vector<std::shared_ptr<const dealii::Triangulation<dim>>> coarse_triangulations;

  /**
   * \brief Collection of DoFhandlers for each multigrid level.
   */
  dealii::MGLevelObject<dealii::DoFHandler<dim>> mg_dof_handlers;

  /**
   * \brief Collection of constraints for each multigrid level.
   */
  dealii::MGLevelObject<dealii::AffineConstraints<float>> level_constraints;

  /**
   * \brief Collection of transfer operators for each multigrid level.
   */
  dealii::MGLevelObject<dealii::MGTwoLevelTransfer<dim, MGVectorType>>
    mg_transfer_operators;

  /**
   * \brief PDE operator for each multigrid level.
   */
  std::unique_ptr<dealii::MGLevelObject<LevelMatrixType>> mg_operators;

  /**
   * \brief Multigrid object for storing all operators.
   */
  std::shared_ptr<dealii::mg::Matrix<MGVectorType>> mg_matrix;

  /**
   * \brief Transfer operator for global coarsening.
   */
  std::shared_ptr<dealii::MGTransferGlobalCoarsening<dim, MGVectorType>> mg_transfer;

  /**
   * \brief Subset of fields that are necessary for the source of the newton update for
   * each multigrid level.
   */
  dealii::MGLevelObject<std::vector<MGVectorType *>> mg_newton_update_src;
};

template <int dim, int degree>
GMGSolver<dim, degree>::GMGSolver(
  const userInputParameters<dim>                       &_user_inputs,
  const variableAttributes                             &_variable_attributes,
  const matrixfreeHandler<dim>                         &_matrix_free_handler,
  const constraintHandler<dim>                         &_constraint_handler,
  const triangulationHandler<dim>                      &_triangulation_handler,
  const prisms::dofHandler<dim>                        &_dof_handler,
  dealii::MGLevelObject<matrixfreeHandler<dim, float>> &_mg_matrix_free_handler,
  solutionHandler<dim>                                 &_solution_handler)
  : linearSolverBase<dim, degree>(_user_inputs,
                                  _variable_attributes,
                                  _matrix_free_handler,
                                  _constraint_handler,
                                  _solution_handler)
  , triangulation_handler(_triangulation_handler)
  , dof_handler(_dof_handler)
  , mg_matrix_free_handler(_mg_matrix_free_handler)
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
  const auto *current_dof_handler = dof_handler.const_dof_handlers.at(this->field_index);

  Assert(this->triangulation_handler.get_triangulation().n_global_levels() > 1,
         dealii::ExcMessage(
           "Multigrid preconditioner requires a multilevel triangulation"));
  if constexpr (dim != 1)
    {
      Assert(this->triangulation_handler.get_triangulation()
               .is_multilevel_hierarchy_constructed(),
             dealii::ExcMessage(
               "parrallel::distributed::Triangulation must be constructed with "
               "multilevel hierarchy"));
    }

  // Create the triangulations for the coarser levels
  coarse_triangulations =
    dealii::MGTransferGlobalCoarseningTools::create_geometric_coarsening_sequence(
      this->triangulation_handler.get_triangulation());

  // Placeholder for p-multigrid
  const auto polynomial_coarsening_sequence =
    dealii::MGTransferGlobalCoarseningTools::create_polynomial_coarsening_sequence(
      current_dof_handler->get_fe().degree,
      dealii::MGTransferGlobalCoarseningTools::PolynomialCoarseningSequenceType::
        decrease_by_one);

  // Create the levels where h is the refinement and p is the polynomial degree.
  std::vector<std::pair<unsigned int, unsigned int>> levels;
  for (unsigned int i = 0; i < coarse_triangulations.size(); ++i)
    {
      levels.emplace_back(i, polynomial_coarsening_sequence.back());
    }

  // Set the maximum and minimum levels for the multigrid based on the triangulation.
  min_level = 0;
  max_level = levels.size() - 1;

  // Init the multilevel operator objects
  mg_operators =
    std::make_unique<dealii::MGLevelObject<LevelMatrixType>>(min_level,
                                                             max_level,
                                                             this->user_inputs,
                                                             this->subset_attributes,
                                                             this->field_index);
  mg_transfer_operators.resize(min_level, max_level);

  // Object for constraints on different levels
  level_constraints.resize(min_level, max_level);
  mg_matrix_free_handler.resize(min_level, max_level, this->user_inputs);

  // Distribute DoFs for each level of the triangulation
  mg_dof_handlers.resize(min_level, max_level);
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      mg_dof_handlers[level].reinit(*coarse_triangulations[level]);

      mg_dof_handlers[level].distribute_dofs(current_dof_handler->get_fe());
    }

  // Apply constraints on each level
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      level_constraints[level].clear();

      const dealii::IndexSet locally_relevant_dofs =
        dealii::DoFTools::extract_locally_relevant_dofs(mg_dof_handlers[level]);

      level_constraints[level].reinit(locally_relevant_dofs);

      dealii::DoFTools::make_hanging_node_constraints(mg_dof_handlers[level],
                                                      level_constraints[level]);

      // TODO: Fix for pinned points
      const auto &boundary_condition =
        this->user_inputs.boundary_parameters.boundary_condition_list.at(
          this->field_index);
      for (const auto &[component, condition] : boundary_condition)
        {
          for (const auto &[boundary_id, boundary_type] :
               condition.boundary_condition_map)
            {
              // Create a mask. This is only applied for vector fields to apply boundary
              // conditions to
              // each component of the vector.
              std::vector<bool> mask(dim, false);
              mask.at(component) = true;

              if (boundary_type == boundaryCondition::type::NATURAL)
                {
                  // Do nothing because they are naturally enforced.
                  continue;
                }
              else if (boundary_type == boundaryCondition::type::DIRICHLET)
                {
                  if (this->user_inputs.var_attributes.at(this->field_index).field_type !=
                      fieldType::VECTOR)
                    {
                      dealii::VectorTools::interpolate_boundary_values(
                        mapping,
                        mg_dof_handlers[level],
                        boundary_id,
                        dealii::Functions::ZeroFunction<dim, float>(1),
                        level_constraints[level]);
                    }
                  else
                    {
                      dealii::VectorTools::interpolate_boundary_values(
                        mapping,
                        mg_dof_handlers[level],
                        boundary_id,
                        dealii::Functions::ZeroFunction<dim, float>(dim),
                        level_constraints[level],
                        mask);
                    }
                }
              else if (boundary_type == boundaryCondition::type::PERIODIC)
                {
                  // Skip boundary ids that are odd since those map to the even faces
                  if (boundary_id % 2 != 0)
                    {
                      continue;
                    }
                  // Create a vector of matched pairs that we fill and enforce upon the
                  // constaints
                  std::vector<dealii::GridTools::PeriodicFacePair<
                    typename dealii::DoFHandler<dim>::cell_iterator>>
                    periodicity_vector;

                  // Determine the direction
                  const auto direction =
                    static_cast<unsigned int>(std::floor(boundary_id / dim));

                  // Collect the matched pairs on the coarsest level of the mesh
                  dealii::GridTools::collect_periodic_faces(mg_dof_handlers[level],
                                                            boundary_id,
                                                            boundary_id + 1,
                                                            direction,
                                                            periodicity_vector);

                  // Set constraints
                  if (this->user_inputs.var_attributes.at(this->field_index).field_type !=
                      fieldType::VECTOR)
                    {
                      dealii::DoFTools::make_periodicity_constraints<dim, dim>(
                        periodicity_vector,
                        level_constraints[level]);
                    }
                  else
                    {
                      dealii::DoFTools::make_periodicity_constraints<dim, dim>(
                        periodicity_vector,
                        level_constraints[level],
                        mask);
                    }
                }
              else if (boundary_type == boundaryCondition::type::NEUMANN)
                {
                  Assert(false,
                         prisms::FeatureNotImplemented("Neumann boundary conditions"));
                }
              else if (boundary_type == boundaryCondition::type::NON_UNIFORM_DIRICHLET)
                {
                  if (this->user_inputs.var_attributes.at(this->field_index).field_type !=
                      fieldType::VECTOR)
                    {
                      dealii::VectorTools::interpolate_boundary_values(
                        mapping,
                        mg_dof_handlers[level],
                        boundary_id,
                        dealii::Functions::ZeroFunction<dim, float>(1),
                        level_constraints[level]);
                    }
                  else
                    {
                      dealii::VectorTools::interpolate_boundary_values(
                        mapping,
                        mg_dof_handlers[level],
                        boundary_id,
                        dealii::Functions::ZeroFunction<dim, float>(dim),
                        level_constraints[level],
                        mask);
                    }
                }
              else if (boundary_type == boundaryCondition::type::NON_UNIFORM_NEUMANN)
                {
                  Assert(false,
                         prisms::FeatureNotImplemented(
                           "Nonuniform neumann boundary conditions"));
                }
            }
        }
      level_constraints[level].close();
    }

  // Setup operator on each level
  mg_newton_update_src.resize(min_level, max_level);
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      // TODO: Fix so mapping is same as rest of the problem. Do the same for the finite
      // element I think.
      // TODO: Fix so that we include all DoF handlers and constraints and select only the
      // ones we need.
      mg_matrix_free_handler[level].reinit(mapping,
                                           mg_dof_handlers[level],
                                           level_constraints[level],
                                           dealii::QGaussLobatto<1>(degree + 1));

      (*mg_operators)[level].initialize(mg_matrix_free_handler[level].get_matrix_free());

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
      // TODO: Check this for all of the vectors. May also be able to remove this
      Assert((*mg_operators)[level]
               .get_matrix_free()
               ->get_vector_partitioner(this->field_index)
               ->is_compatible(*(mg_newton_update_src[level][0]->get_partitioner())),
             dealii::ExcMessage("Incompatabile vector partitioners"));

      // TODO: Allow for src subset that is float instead of double. Not sure how to
      // handle this to be honest. I would rather not store more matrices than I need
      // to.
      (*mg_operators)[level].add_src_solution_subset(mg_newton_update_src[level]);
    }
  mg_matrix = std::make_shared<dealii::mg::Matrix<MGVectorType>>(*mg_operators);

  // Setup transfer operators
  for (unsigned int level = min_level; level < max_level; ++level)
    {
      mg_transfer_operators[level + 1].reinit(mg_dof_handlers[level + 1],
                                              mg_dof_handlers[level],
                                              level_constraints[level + 1],
                                              level_constraints[level]);
    }
  mg_transfer = std::make_shared<dealii::MGTransferGlobalCoarsening<dim, MGVectorType>>(
    mg_transfer_operators);

  this->system_matrix->clear();
  this->system_matrix->initialize(this->matrix_free_handler.get_matrix_free());
  this->update_system_matrix->clear();
  this->update_system_matrix->initialize(this->matrix_free_handler.get_matrix_free());

  this->system_matrix->add_global_to_local_mapping(
    this->residual_global_to_local_solution);
  this->system_matrix->add_src_solution_subset(this->residual_src);

  this->update_system_matrix->add_global_to_local_mapping(
    this->newton_update_global_to_local_solution);
  this->update_system_matrix->add_src_solution_subset(this->newton_update_src);

  // Apply constraints
  this->constraint_handler.get_constraint(this->field_index)
    .distribute(*(this->solution_handler.solution_set.at(
      std::make_pair(this->field_index, dependencyType::NORMAL))));

#ifdef DEBUG
  prisms::conditionalOStreams::pout_summary()
    << "\nMultigrid Setup Information for index " << this->field_index << ":\n"
    << "  Min level: " << min_level << "\n"
    << "  Max level: " << max_level << "\n";
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      prisms::conditionalOStreams::pout_summary()
        << "  Level: " << level << "\n"
        << "    Cells: " << coarse_triangulations[level]->n_global_active_cells() << "\n"
        << "    DoFs: " << mg_dof_handlers[level].n_dofs() << "\n"
        << "    Constrained DoFs: " << level_constraints[level].n_constraints() << "\n";
    }
  prisms::conditionalOStreams::pout_summary()
    << "  MG vertical communication efficiency: "
    << dealii::MGTools::vertical_communication_efficiency(coarse_triangulations) << "\n"
    << "  MG workload imbalance: "
    << dealii::MGTools::workload_imbalance(coarse_triangulations) << "\n\n"
    << std::flush;
#endif
}

template <int dim, int degree>
inline void
GMGSolver<dim, degree>::reinit()
{}

template <int dim, int degree>
inline void
GMGSolver<dim, degree>::solve(const double step_length)
{
  const auto *current_dof_handler = dof_handler.const_dof_handlers.at(this->field_index);
  auto       *solution            = this->solution_handler.solution_set.at(
    std::make_pair(this->field_index, dependencyType::NORMAL));

  // Compute the residual
  this->system_matrix->compute_residual(*this->residual, *solution);
  prisms::conditionalOStreams::pout_summary()
    << "  field: " << this->field_index
    << " Initial residual: " << this->residual->l2_norm() << std::flush;

  // Determine the residual tolerance
  this->compute_solver_tolerance();

  // Update solver controls
  this->solver_control.set_tolerance(this->tolerance);
  dealii::SolverCG<VectorType> cg(this->solver_control);

  // Interpolate the newton update src vector to each multigrid level
  for (unsigned int local_index = 0; local_index < this->newton_update_src.size();
       local_index++)
    {
      // Create a temporary collection of the the dst pointers
      dealii::MGLevelObject<MGVectorType> mg_src_subset(min_level, max_level);
      for (unsigned int level = min_level; level < max_level; ++level)
        {
          mg_src_subset[level] = *mg_newton_update_src[level][local_index];
        }

      // Interpolate
      mg_transfer->interpolate_to_mg(*current_dof_handler,
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
        this->user_inputs.linear_solve_parameters.linear_solve.at(this->field_index)
          .smoothing_range;
      smoother_data[level].degree =
        this->user_inputs.linear_solve_parameters.linear_solve.at(this->field_index)
          .smoother_degree;
      smoother_data[level].eig_cg_n_iterations =
        this->user_inputs.linear_solve_parameters.linear_solve.at(this->field_index)
          .eig_cg_n_iterations;
      (*mg_operators)[level].compute_diagonal(this->field_index);
      smoother_data[level].preconditioner =
        (*mg_operators)[level].get_matrix_diagonal_inverse();
      smoother_data[level].constraints.copy_from(level_constraints[level]);
    }
  mg_smoother.initialize(*mg_operators, smoother_data);

  dealii::MGCoarseGridApplySmoother<MGVectorType> mg_coarse;
  mg_coarse.initialize(mg_smoother);

  // Create multigrid object
  dealii::Multigrid<MGVectorType> mg(*mg_matrix,
                                     mg_coarse,
                                     *mg_transfer,
                                     mg_smoother,
                                     mg_smoother,
                                     min_level,
                                     max_level,
                                     dealii::Multigrid<MGVectorType>::Cycle::v_cycle);

  // Create the preconditioner
  dealii::PreconditionMG<dim,
                         MGVectorType,
                         dealii::MGTransferGlobalCoarsening<dim, MGVectorType>>
    preconditioner(*current_dof_handler, mg, *mg_transfer);

  try
    {
      *this->newton_update = 0.0;
      cg.solve(*(this->update_system_matrix),
               *this->newton_update,
               *this->residual,
               preconditioner);
    }
  catch (...)
    {
      prisms::conditionalOStreams::pout_base()
        << "Warning: linear solver did not converge as per set tolerances.\n";
    }
  this->constraint_handler.get_constraint(this->field_index)
    .set_zero(*this->newton_update);

  prisms::conditionalOStreams::pout_summary()
    << " Final residual: " << this->solver_control.last_value()
    << " Steps: " << this->solver_control.last_step() << "\n"
    << std::flush;

  // Update the solutions
  (*solution).add(step_length, *this->newton_update);
  this->solution_handler.update(fieldSolveType::NONEXPLICIT_LINEAR, this->field_index);

  // Apply constraints
  // This may be redundant with the constraints on the update step.
  this->constraint_handler.get_constraint(this->field_index).distribute(*solution);
}

PRISMS_PF_END_NAMESPACE

#endif