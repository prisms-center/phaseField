#ifndef explicit_postprocess_solver_h
#define explicit_postprocess_solver_h

#include <prismspf/config.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/solvers/explicit_constant_solver.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * Forward declaration for user-implemented PDE class.
 */
template <int dim, int degree, typename number>
class customPDE;

/**
 * \brief This class handles the explicit solves of all postprocessed fields
 */
template <int dim, int degree>
class explicitPostprocessSolver : public explicitBase<dim, degree>
{
public:
  using SystemMatrixType = customPDE<dim, degree, double>;

  /**
   * \brief Constructor.
   */
  explicitPostprocessSolver(const userInputParameters<dim> &_user_inputs,
                            const matrixfreeHandler<dim>   &_matrix_free_handler,
                            const invmHandler<dim, degree> &_invm_handler,
                            const constraintHandler<dim>   &_constraint_handler,
                            const prisms::dofHandler<dim>  &_dof_handler,
                            const dealii::MappingQ1<dim>   &_mapping,
                            solutionHandler<dim>           &_solution_handler);

  /**
   * \brief Destructor.
   */
  ~explicitPostprocessSolver() = default;

  /**
   * \brief Initialize system.
   */
  void
  init() override;

  /**
   * \brief Solve a single update step.
   */
  void
  solve() override;

private:
  /**
   * \brief Mapping from global solution vectors to the local ones
   */
  std::unordered_map<std::pair<unsigned int, dependencyType>, unsigned int, pairHash>
    global_to_local_solution;

  /**
   * \brief Subset of solutions fields that are necessary for explicit solves.
   */
  std::vector<dealii::LinearAlgebra::distributed::Vector<double> *> solution_subset;

  /**
   * \brief Subset of new solutions fields that are necessary for explicit solves.
   */
  std::vector<dealii::LinearAlgebra::distributed::Vector<double> *> new_solution_subset;
};

template <int dim, int degree>
explicitPostprocessSolver<dim, degree>::explicitPostprocessSolver(
  const userInputParameters<dim> &_user_inputs,
  const matrixfreeHandler<dim>   &_matrix_free_handler,
  const invmHandler<dim, degree> &_invm_handler,
  const constraintHandler<dim>   &_constraint_handler,
  const prisms::dofHandler<dim>  &_dof_handler,
  const dealii::MappingQ1<dim>   &_mapping,
  solutionHandler<dim>           &_solution_handler)
  : explicitBase<dim, degree>(_user_inputs,
                              _matrix_free_handler,
                              _invm_handler,
                              _constraint_handler,
                              _dof_handler,
                              _mapping,
                              _solution_handler)
{}

template <int dim, int degree>
inline void
explicitPostprocessSolver<dim, degree>::init()
{
  this->compute_subset_attributes(fieldSolveType::EXPLICIT_POSTPROCESS);

  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  this->compute_shared_dependencies();

  // Create the implementation of customPDE with the subset of variable attributes
  this->system_matrix =
    std::make_unique<SystemMatrixType>(this->user_inputs, this->subset_attributes);

  // Set up the user-implemented equations and create the residual vectors
  this->system_matrix->clear();
  this->system_matrix->initialize(this->matrix_free_handler.get_matrix_free());

  // Create the subset of solution vectors and add the mapping to customPDE
  for (const auto &[index, map] :
       this->subset_attributes.begin()->second.dependency_set_RHS)
    {
      for (const auto &[dependency_type, field_type] : map)
        {
          const auto pair = std::make_pair(index, dependency_type);

          Assert(this->solution_handler.solution_set.find(pair) !=
                   this->solution_handler.solution_set.end(),
                 dealii::ExcMessage("There is no solution vector for the given index = " +
                                    std::to_string(index) +
                                    " and type = " + to_string(dependency_type)));

          Assert(this->solution_handler.new_solution_set.find(index) !=
                   this->solution_handler.new_solution_set.end(),
                 dealii::ExcMessage(
                   "There is no new solution vector for the given index = " +
                   std::to_string(index)));

          solution_subset.push_back(this->solution_handler.solution_set.at(pair));
          new_solution_subset.push_back(
            this->solution_handler.new_solution_set.at(index));
          global_to_local_solution.emplace(pair, solution_subset.size() - 1);
        }
    }

  this->system_matrix->add_global_to_local_mapping(global_to_local_solution);
}

template <int dim, int degree>
inline void
explicitPostprocessSolver<dim, degree>::solve()
{
  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  // Compute the postprocessed fields
  this->system_matrix->compute_postprocess_explicit_update(new_solution_subset,
                                                           solution_subset);

  // Scale the update by the respective (SCALAR/VECTOR) invm. Note that we do this with
  // the original solution set to avoid some messy mapping.
  for (auto [index, vector] : this->solution_handler.new_solution_set)
    {
      if (this->subset_attributes.find(index) != this->subset_attributes.end())
        {
          vector->scale(this->invm_handler.get_invm(index));
        }
    }

  // Update the solutions
  this->solution_handler.update(fieldSolveType::EXPLICIT_POSTPROCESS);
}

PRISMS_PF_END_NAMESPACE

#endif