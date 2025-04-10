#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/explicit_base.h>
#include <prismspf/solvers/explicit_postprocess_solver.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, int degree>
explicitPostprocessSolver<dim, degree>::explicitPostprocessSolver(
  const userInputParameters<dim>                         &_user_inputs,
  const matrixfreeHandler<dim>                           &_matrix_free_handler,
  const invmHandler<dim, degree>                         &_invm_handler,
  const constraintHandler<dim>                           &_constraint_handler,
  const dofHandler<dim>                                  &_dof_handler,
  const dealii::MappingQ1<dim>                           &_mapping,
  solutionHandler<dim>                                   &_solution_handler,
  std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator)
  : explicitBase<dim, degree>(_user_inputs,
                              _matrix_free_handler,
                              _invm_handler,
                              _constraint_handler,
                              _dof_handler,
                              _mapping,
                              _solution_handler,
                              std::move(_pde_operator))
{}

template <int dim, int degree>
void
explicitPostprocessSolver<dim, degree>::init()
{
  this->compute_subset_attributes(fieldSolveType::EXPLICIT_POSTPROCESS);

  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  this->compute_shared_dependencies();

  // Create the implementation of matrixFreeOperator with the subset of variable
  // attributes
  this->system_matrix =
    std::make_unique<SystemMatrixType>(this->subset_attributes, this->pde_operator);

  // Set up the user-implemented equations and create the residual vectors
  this->system_matrix->clear();
  this->system_matrix->initialize(this->matrix_free_handler->get_matrix_free());

  // Create the subset of solution vectors and add the mapping to matrixFreeOperator
  for (const auto &[index, map] :
       this->subset_attributes.begin()->second.dependency_set_RHS)
    {
      for (const auto &[dependency_type, field_type] : map)
        {
          const auto pair = std::make_pair(index, dependency_type);

          solution_subset.push_back(
            this->solution_handler->get_solution_vector(index, dependency_type));
          new_solution_subset.push_back(
            this->solution_handler->get_new_solution_vector(index));
          global_to_local_solution.emplace(pair, solution_subset.size() - 1);
        }
    }

  this->system_matrix->add_global_to_local_mapping(global_to_local_solution);
}

template <int dim, int degree>
void
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
  for (auto [index, vector] : this->solution_handler->get_new_solution_vector())
    {
      if (this->subset_attributes.find(index) != this->subset_attributes.end())
        {
          vector->scale(this->invm_handler->get_invm(index));
        }
    }

  // Update the solutions
  this->solution_handler->update(fieldSolveType::EXPLICIT_POSTPROCESS);
}

INSTANTIATE_BI_TEMPLATE(explicitPostprocessSolver)

PRISMS_PF_END_NAMESPACE