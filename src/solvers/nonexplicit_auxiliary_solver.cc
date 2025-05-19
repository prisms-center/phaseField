// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/mg_level_object.h>
#include <deal.II/fe/mapping_q1.h>

#include <prismspf/core/constraint_handler.h>
#include <prismspf/core/dof_handler.h>
#include <prismspf/core/invm_handler.h>
#include <prismspf/core/matrix_free_handler.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/solution_handler.h>
#include <prismspf/core/triangulation_handler.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/nonexplicit_auxiliary_solver.h>
#include <prismspf/solvers/nonexplicit_base.h>

#include <prismspf/config.h>

#include <map>
#include <memory>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree>
nonexplicitAuxiliarySolver<dim, degree>::nonexplicitAuxiliarySolver(
  const userInputParameters<dim>                         &_user_inputs,
  const matrixfreeHandler<dim>                           &_matrix_free_handler,
  const triangulationHandler<dim>                        &_triangulation_handler,
  const invmHandler<dim, degree>                         &_invm_handler,
  const constraintHandler<dim, degree>                   &_constraint_handler,
  const dofHandler<dim>                                  &_dof_handler,
  const dealii::MappingQ1<dim>                           &_mapping,
  dealii::MGLevelObject<matrixfreeHandler<dim, float>>   &_mg_matrix_free_handler,
  solutionHandler<dim>                                   &_solution_handler,
  std::shared_ptr<const PDEOperator<dim, degree, double>> _pde_operator)
  : nonexplicitBase<dim, degree>(_user_inputs,
                                 _matrix_free_handler,
                                 _triangulation_handler,
                                 _invm_handler,
                                 _constraint_handler,
                                 _dof_handler,
                                 _mapping,
                                 _mg_matrix_free_handler,
                                 _solution_handler,
                                 std::move(_pde_operator))
{}

template <unsigned int dim, unsigned int degree>
inline void
nonexplicitAuxiliarySolver<dim, degree>::init()
{
  this->compute_subset_attributes(fieldSolveType::NONEXPLICIT_AUXILIARY);

  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }

  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      // Creating temporary map to match types
      std::map<unsigned int, variableAttributes> temp;
      temp.emplace(index, variable);
      subset_attributes_list.push_back(temp);

      // Create the implementation of matrixFreeOperator with the subset of variable
      // attributes
      this->get_system_matrix()[index] =
        std::make_unique<SystemMatrixType>(this->get_subset_attributes(),
                                           this->get_pde_operator(),
                                           index);

      // Set up the user-implemented equations and create the residual vectors
      this->get_system_matrix()[index]->clear();
      this->get_system_matrix()[index]->initialize(
        this->get_matrix_free_handler().get_matrix_free());

      // Create the subset of solution vectors and add the mapping to matrixFreeOperator
      new_solution_subset[index].push_back(
        this->get_solution_handler().get_new_solution_vector(index));
      solution_subset[index].push_back(
        this->get_solution_handler().get_solution_vector(index, dependencyType::NORMAL));
      global_to_local_solution[index].emplace(std::make_pair(index,
                                                             dependencyType::NORMAL),
                                              0);
      for (const auto &[variable_index, map] :
           this->get_subset_attributes().begin()->second.get_dependency_set_rhs())
        {
          for (const auto &[dependency_type, field_type] : map)
            {
              const auto pair = std::make_pair(variable_index, dependency_type);

              solution_subset[index].push_back(
                this->get_solution_handler().get_solution_vector(variable_index,
                                                                 dependency_type));
              global_to_local_solution[index].emplace(pair,
                                                      solution_subset.at(index).size() -
                                                        1);
            }
        }
      this->get_system_matrix()[index]->add_global_to_local_mapping(
        global_to_local_solution.at(index));
    }
}

template <unsigned int dim, unsigned int degree>
inline void
nonexplicitAuxiliarySolver<dim, degree>::solve()
{
  // If the subset attribute is empty return early
  if (this->get_subset_attributes().empty())
    {
      return;
    }

  for (const auto &[index, variable] : this->get_subset_attributes())
    {
      // Compute the update
      this->get_system_matrix()[index]->compute_nonexplicit_auxiliary_update(
        new_solution_subset.at(index),
        solution_subset.at(index));

      // Scale the update by the respective (SCALAR/VECTOR) invm.
      new_solution_subset.at(index).at(0)->scale(
        this->get_invm_handler().get_invm(index));

      // Update the solutions
      this->get_solution_handler().update(fieldSolveType::NONEXPLICIT_AUXILIARY, index);

      // Apply constraints
      this->get_constraint_handler().get_constraint(index).distribute(*(
        this->get_solution_handler().get_solution_vector(index, dependencyType::NORMAL)));
    }
}

INSTANTIATE_BI_TEMPLATE(nonexplicitAuxiliarySolver)

PRISMS_PF_END_NAMESPACE