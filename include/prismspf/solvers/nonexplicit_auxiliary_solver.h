// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#ifndef nonexplicit_auxiliary_solver_h
#define nonexplicit_auxiliary_solver_h

#include <prismspf/config.h>
#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/solvers/nonexplicit_base.h>
#include <prismspf/user_inputs/user_input_parameters.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * Forward declaration for user-implemented PDE class.
 */
template <int dim, int degree, typename number>
class customPDE;

/**
 * \brief This class handles all auxiliary solves.
 */
template <int dim, int degree>
class nonexplicitAuxiliarySolver : public nonexplicitBase<dim, degree>
{
public:
  using SystemMatrixType = customPDE<dim, degree, double>;
  using VectorType       = dealii::LinearAlgebra::distributed::Vector<double>;

  /**
   * \brief Constructor.
   */
  nonexplicitAuxiliarySolver(
    const userInputParameters<dim>                       &_user_inputs,
    const matrixfreeHandler<dim>                         &_matrix_free_handler,
    const triangulationHandler<dim>                      &_triangulation_handler,
    const invmHandler<dim, degree>                       &_invm_handler,
    const constraintHandler<dim>                         &_constraint_handler,
    const prisms::dofHandler<dim>                        &_dof_handler,
    const dealii::MappingQ1<dim>                         &_mapping,
    dealii::MGLevelObject<matrixfreeHandler<dim, float>> &_mg_matrix_free_handler,
    solutionHandler<dim>                                 &_solution_handler);

  /**
   * \brief Destructor.
   */
  ~nonexplicitAuxiliarySolver() = default;

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
  std::map<
    unsigned int,
    std::unordered_map<std::pair<unsigned int, dependencyType>, unsigned int, pairHash>>
    global_to_local_solution;

  /**
   * \brief Subset of solutions fields that are necessary for explicit solves.
   */
  std::map<unsigned int, std::vector<VectorType *>> solution_subset;

  /**
   * \brief Subset of new solutions fields that are necessary for explicit solves.
   */
  std::map<unsigned int, std::vector<VectorType *>> new_solution_subset;

  /**
   * \brief List of subset attributes.
   */
  std::vector<std::map<unsigned int, variableAttributes>> subset_attributes_list;
};

template <int dim, int degree>
nonexplicitAuxiliarySolver<dim, degree>::nonexplicitAuxiliarySolver(
  const userInputParameters<dim>                       &_user_inputs,
  const matrixfreeHandler<dim>                         &_matrix_free_handler,
  const triangulationHandler<dim>                      &_triangulation_handler,
  const invmHandler<dim, degree>                       &_invm_handler,
  const constraintHandler<dim>                         &_constraint_handler,
  const prisms::dofHandler<dim>                        &_dof_handler,
  const dealii::MappingQ1<dim>                         &_mapping,
  dealii::MGLevelObject<matrixfreeHandler<dim, float>> &_mg_matrix_free_handler,
  solutionHandler<dim>                                 &_solution_handler)
  : nonexplicitBase<dim, degree>(_user_inputs,
                                 _matrix_free_handler,
                                 _triangulation_handler,
                                 _invm_handler,
                                 _constraint_handler,
                                 _dof_handler,
                                 _mapping,
                                 _mg_matrix_free_handler,
                                 _solution_handler)
{}

template <int dim, int degree>
inline void
nonexplicitAuxiliarySolver<dim, degree>::init()
{
  this->compute_subset_attributes(fieldSolveType::NONEXPLICIT_AUXILIARY);

  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  for (const auto &[index, variable] : this->subset_attributes)
    {
      // Creating temporary map to match types
      std::map<unsigned int, variableAttributes> temp;
      temp.emplace(index, variable);
      subset_attributes_list.push_back(temp);

      // Create the implementation of customPDE with the subset of variable attributes
      this->system_matrix[index] =
        std::make_unique<SystemMatrixType>(this->user_inputs,
                                           index,
                                           subset_attributes_list.back());

      // Set up the user-implemented equations and create the residual vectors
      this->system_matrix.at(index)->clear();
      this->system_matrix.at(index)->initialize(
        this->matrix_free_handler.get_matrix_free());

      // Create the subset of solution vectors and add the mapping to customPDE
      new_solution_subset[index].push_back(
        this->solution_handler.new_solution_set.at(index));
      solution_subset[index].push_back(this->solution_handler.solution_set.at(
        std::make_pair(index, dependencyType::NORMAL)));
      global_to_local_solution[index].emplace(std::make_pair(index,
                                                             dependencyType::NORMAL),
                                              0);
      for (const auto &[variable_index, map] :
           subset_attributes_list.back().begin()->second.dependency_set_RHS)
        {
          for (const auto &[dependency_type, field_type] : map)
            {
              const auto pair = std::make_pair(variable_index, dependency_type);

              Assert(this->solution_handler.solution_set.find(pair) !=
                       this->solution_handler.solution_set.end(),
                     dealii::ExcMessage(
                       "There is no solution vector for the given index = " +
                       std::to_string(variable_index) +
                       " and type = " + to_string(dependency_type)));

              Assert(this->solution_handler.new_solution_set.find(variable_index) !=
                       this->solution_handler.new_solution_set.end(),
                     dealii::ExcMessage(
                       "There is no new solution vector for the given index = " +
                       std::to_string(variable_index)));

              solution_subset[index].push_back(
                this->solution_handler.solution_set.at(pair));
              global_to_local_solution[index].emplace(pair,
                                                      solution_subset.at(index).size() -
                                                        1);
            }
        }
      this->system_matrix.at(index)->add_global_to_local_mapping(
        global_to_local_solution.at(index));
    }
}

template <int dim, int degree>
inline void
nonexplicitAuxiliarySolver<dim, degree>::solve()
{
  // If the subset attribute is empty return early
  if (this->subset_attributes.empty())
    {
      return;
    }

  for (const auto &[index, variable] : this->subset_attributes)
    {
      // Compute the update
      this->system_matrix.at(index)->compute_nonexplicit_auxiliary_update(
        new_solution_subset.at(index),
        solution_subset.at(index));

      // Scale the update by the respective (SCALAR/VECTOR) invm.
      new_solution_subset.at(index).at(0)->scale(this->invm_handler.get_invm(index));

      // Update the solutions
      this->solution_handler.update(fieldSolveType::NONEXPLICIT_AUXILIARY, index);

      // Apply constraints
      this->constraint_handler.get_constraint(index).distribute(
        *(this->solution_handler.solution_set.at(
          std::make_pair(index, dependencyType::NORMAL))));
    }
}

PRISMS_PF_END_NAMESPACE

#endif