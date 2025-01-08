#ifndef gmg_solver_h
#define gmg_solver_h

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/operators.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>
#include <deal.II/multigrid/multigrid.h>

#include <core/solvers/base_solver.h>

/**
 * Forward declaration for user-implemented PDE class.
 */
template <int dim, int degree, typename number>
class customPDE;

/**
 * \brief Class that handles the assembly and solving of a field with a GMG preconditioner
 */
template <int dim, int degree>
class GMGSolver : public baseSolver<dim, degree>
{
public:
  using SystemMatrixType = customPDE<dim, degree, double>;

  using LevelMatrixType = customPDE<dim, degree, float>;

  /**
   * \brief Constructor.
   */
  GMGSolver(const userInputParameters<dim> &_user_inputs, const uint &_field_index);

  /**
   * \brief Destructor.
   */
  ~GMGSolver() = default;

  /**
   * \brief Initialize the system.
   */
  void
  init_system(const triangulationHandler<dim>                    &triangulation_handler,
              dealii::LinearAlgebra::distributed::Vector<double> &solution,
              dealii::DoFHandler<dim>                            &dof_handler) override;

  /**
   * \brief Reinitialize the system.
   */
  void
  reinit_system() override;

  /**
   * \brief Solve the system Ax=b.
   */
  void
  solve(const triangulationHandler<dim>                    &triangulation_handler,
        dealii::LinearAlgebra::distributed::Vector<double> &solution,
        const dealii::DoFHandler<dim>                      &dof_handler) override;

private:
  /**
   * \brief PDE operator for each multigrid level.
   */
  dealii::MGLevelObject<LevelMatrixType> mg_matrices;

  /**
   * \brief Boundary constraints for multigrid levels.
   */
  dealii::MGConstrainedDoFs mg_constrained_dofs;

  /**
   * \brief Set of boundary ids that are dirichlet constraints
   */
  std::set<dealii::types::boundary_id> dirichlet_boundary_ids;
};

template <int dim, int degree>
GMGSolver<dim, degree>::GMGSolver(const userInputParameters<dim> &_user_inputs,
                                  const uint                     &_field_index)
  : baseSolver<dim, degree>(_user_inputs, _field_index)
  , mg_matrices(dealii::MGLevelObject<LevelMatrixType>(0, 0, _user_inputs))
{
  // Get the dirichlet boundary ids for multigrid constriants
  for (const auto &[component, boundary_condition] :
       this->user_inputs.boundary_parameters.boundary_condition_list.at(
         this->field_index))
    {
      Assert(component == 0, FeatureNotImplemented("Vector multigrid"));

      for (const auto &[boundary_id, boundary_type] :
           boundary_condition.boundary_condition_map)
        {
          if (boundary_type == boundaryType::DIRICHLET ||
              boundary_type == boundaryType::NON_UNIFORM_DIRICHLET)
            {
              dirichlet_boundary_ids.insert(boundary_id);
            }
        }
    }
}

template <int dim, int degree>
inline void
GMGSolver<dim, degree>::init_system(
  const triangulationHandler<dim>                    &triangulation_handler,
  dealii::LinearAlgebra::distributed::Vector<double> &solution,
  dealii::DoFHandler<dim>                            &dof_handler)
{
  this->system_matrix.clear();
  mg_matrices.clear_elements();

  dof_handler.reinit(triangulation_handler.get_triangulation());
  dof_handler.distribute_dofs(this->fe);
  dof_handler.distribute_mg_dofs();

  this->constraint_handler.make_constraints(this->mapping, dof_handler);

  typename dealii::MatrixFree<dim, double>::AdditionalData additional_data;
  additional_data.tasks_parallel_scheme =
    dealii::MatrixFree<dim, double>::AdditionalData::none;

  // This should be done according to the residual flags to prevent excess data being
  // evaluated for the shape functions.
  additional_data.mapping_update_flags =
    (update_values | update_gradients | update_JxW_values | update_quadrature_points);
  std::shared_ptr<dealii::MatrixFree<dim, double>> system_mf_storage(
    new dealii::MatrixFree<dim, double>());
  system_mf_storage->reinit(this->mapping,
                            dof_handler,
                            this->constraint_handler.get_constraints(),
                            dealii::QGauss<1>(this->fe.degree + 1),
                            additional_data);
  this->system_matrix.initialize(system_mf_storage);

  this->system_matrix.initialize_dof_vector(this->change_in_solution);
  this->system_matrix.initialize_dof_vector(solution);
  this->system_matrix.initialize_dof_vector(this->residual);

  // Apply boundary conditions to the solution vector as an initial guess
  this->constraint_handler.get_constraints().distribute(solution);

  const uint nlevels = triangulation_handler.get_n_global_levels();
  mg_matrices.resize(0, nlevels - 1, this->user_inputs);

  mg_constrained_dofs.initialize(dof_handler);
  mg_constrained_dofs.make_zero_boundary_constraints(dof_handler, dirichlet_boundary_ids);

  for (uint level = 0; level < nlevels; ++level)
    {
      // Setup the constraint set for the multigrid level
      dealii::AffineConstraints<double> level_constraints(
        dof_handler.locally_owned_mg_dofs(level),
        dealii::DoFTools::extract_locally_relevant_level_dofs(dof_handler, level));

      // Loop through boundary indices and apply constraints
      for (const dealii::types::global_dof_index dof_index :
           mg_constrained_dofs.get_boundary_indices(level))
        {
          level_constraints.constrain_dof_to_zero(dof_index);
        }
      level_constraints.close();

      typename dealii::MatrixFree<dim, float>::AdditionalData additional_data;
      additional_data.tasks_parallel_scheme =
        dealii::MatrixFree<dim, float>::AdditionalData::none;

      // This should be done according to the residual flags to prevent excess data being
      // evaluated for the shape functions.
      additional_data.mapping_update_flags =
        (update_values | update_gradients | update_JxW_values | update_quadrature_points);
      additional_data.mg_level = level;
      std::shared_ptr<dealii::MatrixFree<dim, float>> mg_mf_storage_level =
        std::make_shared<dealii::MatrixFree<dim, float>>();
      mg_mf_storage_level->reinit(this->mapping,
                                  dof_handler,
                                  level_constraints,
                                  dealii::QGauss<1>(this->fe.degree + 1),
                                  additional_data);

      mg_matrices[level].initialize(mg_mf_storage_level, mg_constrained_dofs, level);
    }

  conditionalOStreams::pout_base
    << "Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;
  conditionalOStreams::pout_summary()
    << "Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;
}

template <int dim, int degree>
inline void
GMGSolver<dim, degree>::reinit_system()
{}

template <int dim, int degree>
inline void
GMGSolver<dim, degree>::solve(
  const triangulationHandler<dim>                    &triangulation_handler,
  dealii::LinearAlgebra::distributed::Vector<double> &solution,
  const dealii::DoFHandler<dim>                      &dof_handler)
{
  dealii::MGTransferMatrixFree<dim, float> mg_transfer(mg_constrained_dofs);
  mg_transfer.build(dof_handler);

  using SmootherType =
    dealii::PreconditionChebyshev<LevelMatrixType,
                                  dealii::LinearAlgebra::distributed::Vector<float>>;
  dealii::mg::SmootherRelaxation<SmootherType,
                                 dealii::LinearAlgebra::distributed::Vector<float>>
                                                               mg_smoother;
  dealii::MGLevelObject<typename SmootherType::AdditionalData> smoother_data;
  smoother_data.resize(0, triangulation_handler.get_n_global_levels() - 1);
  for (uint level = 0; level < triangulation_handler.get_n_global_levels(); ++level)
    {
      if (level > 0)
        {
          smoother_data[level].smoothing_range =
            this->user_inputs.linear_solve_parameters.linear_solve.at(this->field_index)
              .smoothing_range;
          smoother_data[level].degree =
            this->user_inputs.linear_solve_parameters.linear_solve.at(this->field_index)
              .smoother_iterations;
          smoother_data[level].eig_cg_n_iterations =
            this->user_inputs.linear_solve_parameters.linear_solve.at(this->field_index)
              .eig_cg_n_iterations;
        }
      else
        {
          smoother_data[0].smoothing_range     = 1e-3;
          smoother_data[0].degree              = dealii::numbers::invalid_unsigned_int;
          smoother_data[0].eig_cg_n_iterations = mg_matrices[0].m();
        }
      mg_matrices[level].compute_diagonal();
      smoother_data[level].preconditioner =
        mg_matrices[level].get_matrix_diagonal_inverse();
    }
  mg_smoother.initialize(mg_matrices, smoother_data);

  dealii::MGCoarseGridApplySmoother<dealii::LinearAlgebra::distributed::Vector<float>>
    mg_coarse;
  mg_coarse.initialize(mg_smoother);

  dealii::mg::Matrix<dealii::LinearAlgebra::distributed::Vector<float>> mg_matrix(
    mg_matrices);

  dealii::MGLevelObject<dealii::MatrixFreeOperators::MGInterfaceOperator<LevelMatrixType>>
    mg_interface_matrices;
  mg_interface_matrices.resize(0, triangulation_handler.get_n_global_levels() - 1);
  for (uint level = 0; level < triangulation_handler.get_n_global_levels(); ++level)
    {
      mg_interface_matrices[level].initialize(mg_matrices[level]);
    }
  dealii::mg::Matrix<dealii::LinearAlgebra::distributed::Vector<float>> mg_interface(
    mg_interface_matrices);

  Multigrid<LinearAlgebra::distributed::Vector<float>> mg(mg_matrix,
                                                          mg_coarse,
                                                          mg_transfer,
                                                          mg_smoother,
                                                          mg_smoother);
  mg.set_edge_matrices(mg_interface, mg_interface);

  dealii::PreconditionMG<dim,
                         dealii::LinearAlgebra::distributed::Vector<float>,
                         dealii::MGTransferMatrixFree<dim, float>>
    preconditioner(dof_handler, mg, mg_transfer);

  // Compute the residual of the current solution (Ax-b)
  this->system_matrix.compute_residual(this->residual, solution);

  // Determine the residual tolerance
  this->compute_solver_tolerance();

  // Update solver controls
  this->solver_control.set_tolerance(this->tolerance);
  dealii::SolverCG<dealii::LinearAlgebra::distributed::Vector<double>> cg(
    this->solver_control);

  try
    {
      this->change_in_solution = 0.0;
      cg.solve(this->system_matrix,
               this->change_in_solution,
               this->residual,
               preconditioner);
      this->constraint_handler.get_constraints().set_zero(this->change_in_solution);
    }
  catch (...)
    {
      conditionalOStreams::pout_base
        << "Warning: linear solver did not converge as per set tolerances.\n";
    }
  solution += this->change_in_solution;
  this->constraint_handler.get_constraints().distribute(solution);

  conditionalOStreams::pout_summary()
    << this->solver_control.last_step() << " iterations\n";
}

#endif