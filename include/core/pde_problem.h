#ifndef pde_problem_h
#define pde_problem_h

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <core/boundary_conditions/constraint_handler.h>
#include <core/conditional_ostreams.h>
#include <core/invm_handler.h>
#include <core/matrix_free_handler.h>
#include <core/solution_output.h>
#include <core/triangulation_handler.h>
#include <core/type_enums.h>
#include <core/user_inputs/user_input_parameters.h>
#include <core/variable_attributes.h>
#include <map>
#include <memory>
#include <vector>

/**
 * \brief This is the main class that handles the construction and solving of
 * user-specified PDEs.
 */
template <int dim, int degree>
class PDEProblem
{
public:
  /**
   * \brief Constructor.
   */
  PDEProblem(const userInputParameters<dim> &_user_inputs);

  /**
   * \brief Destructor.
   */
  ~PDEProblem();

  /**
   * \brief Run initialization and solving steps of the given problem.
   */
  void
  run();

private:
  /**
   * \brief Main time-stepping loop that calls solve_increment, reinit_system,
   * output_results, etc...
   */
  void
  solve();

  /**
   * \brief Solve a single increment of the given PDEs.
   */
  void
  solve_increment();

  /**
   * \brief Initialize the system.
   */
  void
  init_system();

  /**
   * \brief Reinitialize the system.
   */
  void
  reinit_system();

  /**
   * \brief User-inputs.
   */
  const userInputParameters<dim> &user_inputs;

  /**
   * \brief Triangulation handler.
   */
  std::shared_ptr<triangulationHandler<dim>> triangulation_handler;

  /**
   * \brief Constraint handler.
   */
  std::shared_ptr<constraintHandler<dim>> constraint_handler;

  /**
   * \brief Matrix-free object handler for non-multigrid data.
   */
  std::shared_ptr<matrixfreeHandler<dim>> matrix_free_handler;

  /**
   * \brief Matrix-free object handler for multigrid data.
   */
  std::shared_ptr<matrixfreeHandler<dim, float>> multigrid_matrix_free_handler;

  /**
   * \brief invm handler.
   */
  std::shared_ptr<invmHandler<dim, degree>> invm_handler;

  /**
   * \brief Collection of the triangulation DoFs. The number of DoFHandlers should be
   * equal to or less than the number of fields. Technically, there's a small optimization
   * we can use when multiple fields have the same constraints and quadrature rule,
   * allowing us to share the same DoFHandler. An example of this might be grain growth.
   */
  std::vector<dealii::DoFHandler<dim> *> dof_handlers;

  /**
   * \brief Const copy of the dof_handlers;
   */
  std::vector<const dealii::DoFHandler<dim> *> const_dof_handlers;

  /**
   * \brief Collection of finite element systems. This is just a collection of two
   * FESystem's: one for scalar fields and one for vector fields. For now they both use
   * FE_Q finite elements.
   */
  std::map<fieldType, dealii::FESystem<dim>> fe_system;

  /**
   * \brief Mappings to and from reference cell.
   */
  dealii::MappingQ1<dim> mapping;

  /**
   * \brief The collection of solution vector at the current timestep. This includes
   * current values and old values.
   */
  std::unordered_map<std::pair<uint, dependencyType>,
                     dealii::LinearAlgebra::distributed::Vector<double> *,
                     pairHash>
    solution_set;

  /**
   * \brief Number of unique fields with different boundary conditions and fieldType
   * (SCALAR/VECTOR).
   */
  uint n_unique_fields = 0;
};

template <int dim, int degree>
PDEProblem<dim, degree>::PDEProblem(const userInputParameters<dim> &_user_inputs)
  : user_inputs(_user_inputs)
  , triangulation_handler(std::make_shared<triangulationHandler<dim>>(_user_inputs))
  , constraint_handler(std::make_shared<constraintHandler<dim>>(_user_inputs))
  , matrix_free_handler(std::make_shared<matrixfreeHandler<dim>>(_user_inputs))
  , multigrid_matrix_free_handler(
      std::make_shared<matrixfreeHandler<dim, float>>(_user_inputs))
  , invm_handler(std::make_shared<invmHandler<dim, degree>>(_user_inputs.var_attributes))
{
  dof_handlers.resize(user_inputs.var_attributes.size());
  for (const auto &[index, variable] : user_inputs.var_attributes)
    {
      dof_handlers.at(index) = new dealii::DoFHandler<dim>();
    }
  for (auto &dof_handler : dof_handlers)
    {
      const_dof_handlers.push_back(dof_handler);
    }
}

template <int dim, int degree>
PDEProblem<dim, degree>::~PDEProblem()
{
  for (auto dof_handler : dof_handlers)
    {
      delete dof_handler;
    }
  dof_handlers.clear();
  for (const auto &[pair, solution] : solution_set)
    {
      delete solution;
    }
  solution_set.clear();
}

template <int dim, int degree>
void
PDEProblem<dim, degree>::init_system()
{
  // Create the SCALAR/VECTOR FESystem's, if applicable
  for (const auto &[index, variable] : user_inputs.var_attributes)
    {
      if (variable.field_type == fieldType::SCALAR &&
          fe_system.find(fieldType::SCALAR) == fe_system.end())
        {
          fe_system.emplace(fieldType::SCALAR,
                            dealii::FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(degree + 1)),
                                                  1));
        }
      else if (variable.field_type == fieldType::VECTOR &&
               fe_system.find(fieldType::VECTOR) == fe_system.end())
        {
          fe_system.emplace(fieldType::VECTOR,
                            dealii::FESystem<dim>(FE_Q<dim>(QGaussLobatto<1>(degree + 1)),
                                                  dim));
        }
    }

  // Create the mesh
  triangulation_handler->generate_mesh();

  // Create the dof handlers.
  // TODO: Fix so we can make unique instances of dof handlers.
  uint n_dofs = 0;
  for (const auto &[index, variable] : user_inputs.var_attributes)
    {
      dof_handlers.at(index)->reinit(triangulation_handler->get_triangulation());
      dof_handlers.at(index)->distribute_dofs(fe_system.at(variable.field_type));

      // If we have GMG enabled distribute those constraints as well to that field.
      if (user_inputs.linear_solve_parameters.linear_solve.find(index) !=
          user_inputs.linear_solve_parameters.linear_solve.end())
        {
          if (user_inputs.linear_solve_parameters.linear_solve.at(index).preconditioner ==
              preconditionerType::GMG)
            {
              dof_handlers.at(index)->distribute_mg_dofs();
            }
        }
      n_dofs += dof_handlers.at(index)->n_dofs();
    }
  conditionalOStreams::pout_base << "Number of degrees of freedom: " << n_dofs << "\n";
  conditionalOStreams::pout_summary()
    << "Number of degrees of freedom: " << n_dofs << "\n";

  // Create the constraints
  constraint_handler->make_constraints(mapping, dof_handlers);

  // Reinit the matrix-free objects
  matrix_free_handler->reinit(mapping,
                              const_dof_handlers,
                              constraint_handler->get_constraints(),
                              QGaussLobatto<1>(degree + 1));

  // Create the solution set
  for (const auto &[index, variable] : user_inputs.var_attributes)
    {
      // Add the current variable if it doesn't already exist
      if (solution_set.find(std::make_pair(index, dependencyType::NORMAL)) ==
          solution_set.end())
        {
          solution_set[std::make_pair(index, dependencyType::NORMAL)] =
            new dealii::LinearAlgebra::distributed::Vector<double>();
        }

      // Add dependencies if they don't exist
      for (const auto &[pair, flags] : variable.eval_flag_set_RHS)
        {
          if (pair.second == dependencyType::CHANGE)
            {
              continue;
            }
          if (solution_set.find(pair) == solution_set.end())
            {
              solution_set[pair] =
                new dealii::LinearAlgebra::distributed::Vector<double>();
            }
        }
      for (const auto &[pair, flags] : variable.eval_flag_set_LHS)
        {
          if (pair.second == dependencyType::CHANGE)
            {
              continue;
            }
          if (solution_set.find(pair) == solution_set.end())
            {
              solution_set[pair] =
                new dealii::LinearAlgebra::distributed::Vector<double>();
            }
        }
    }

  // Initialize the invm and compute it
  invm_handler->initialize(matrix_free_handler->get_matrix_free());
  invm_handler->compute_invm();

  solutionOutput<dim> output(invm_handler->get_invm(0),
                             *(dof_handlers[0]),
                             degree,
                             "invm",
                             0);

  // Initialize the solver types
}

template <int dim, int degree>
void
PDEProblem<dim, degree>::solve_increment()
{}

template <int dim, int degree>
void
PDEProblem<dim, degree>::solve()
{
  init_system();
  // solve_increment();

  // solution.update_ghost_values();
  // solutionOutput<dim> output(solution, dof_handler, degree);
}

template <int dim, int degree>
void
PDEProblem<dim, degree>::run()
{
  // Print information regarding the vectorized array lanes to summary.log
  const uint n_vect_doubles = dealii::VectorizedArray<double>::size();
  const uint n_vect_bits    = 8 * sizeof(double) * n_vect_doubles;

  conditionalOStreams::pout_summary()
    << "================================================\n"
    << "\tSolve\n"
    << "================================================\n"
    << "Vectorization over " << n_vect_doubles << " doubles = " << n_vect_bits
    << " bits (" << Utilities::System::get_current_vectorization_level() << ')' << '\n';

  solve();
}

#endif