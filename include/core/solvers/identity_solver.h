#ifndef identity_solver_h
#define identity_solver_h

#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <core/solvers/base_solver.h>

/**
 * Forward declaration for user-implemented PDE class.
 */
template <int dim, int degree, typename number>
class customPDE;

/**
 * \brief Class that handles the assembly and solving of a field with the identity
 * preconditioner (no preconditioner)
 */
template <int dim, int degree>
class identitySolver : public baseSolver<dim, degree>
{
public:
  using SystemMatrixType = customPDE<dim, degree, double>;

  /**
   * \brief Constructor.
   */
  identitySolver(const userInputParameters<dim> &_user_inputs, const uint &_field_index);

  /**
   * \brief Destructor.
   */
  ~identitySolver() = default;

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
};

template <int dim, int degree>
identitySolver<dim, degree>::identitySolver(const userInputParameters<dim> &_user_inputs,
                                            const uint                     &_field_index)
  : baseSolver<dim, degree>(_user_inputs, _field_index)
{}

template <int dim, int degree>
inline void
identitySolver<dim, degree>::init_system(
  const triangulationHandler<dim>                    &triangulation_handler,
  dealii::LinearAlgebra::distributed::Vector<double> &solution,
  dealii::DoFHandler<dim>                            &dof_handler)
{
  this->system_matrix.clear();

  dof_handler.reinit(triangulation_handler.get_triangulation());
  dof_handler.distribute_dofs(this->fe);

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

  conditionalOStreams::pout_base
    << "Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;
  conditionalOStreams::pout_summary()
    << "Number of degrees of freedom: " << dof_handler.n_dofs() << std::endl;
}

template <int dim, int degree>
inline void
identitySolver<dim, degree>::reinit_system()
{}

template <int dim, int degree>
inline void
identitySolver<dim, degree>::solve(
  [[maybe_unused]] const triangulationHandler<dim>   &triangulation_handler,
  dealii::LinearAlgebra::distributed::Vector<double> &solution,
  [[maybe_unused]] const dealii::DoFHandler<dim>     &dof_handler)
{
  // Compute the residual of the current solution (Ax-b)
  this->system_matrix.compute_residual(this->residual, solution);

  // Determine the residual tolerance
  this->compute_solver_tolerance();

  // Setup solver controls
  this->solver_control.set_tolerance(this->tolerance);
  dealii::SolverCG<dealii::LinearAlgebra::distributed::Vector<double>> cg(
    this->solver_control);

  try
    {
      this->change_in_solution = 0.0;
      cg.solve(this->system_matrix,
               this->change_in_solution,
               this->residual,
               PreconditionIdentity());
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