// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/types.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/exceptions.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/solvers/mf_operator.h>

#include <prismspf/utilities/element_volume.h>

#include <prismspf/config.h>

#include "prismspf/core/group_solution_handler.h"
#include "prismspf/core/solve_group.h"

#include <map>
#include <memory>
#include <utility>
#include <vector>

#if DEAL_II_VERSION_MAJOR >= 9 && DEAL_II_VERSION_MINOR >= 7
#  include <deal.II/base/enable_observer_pointer.h>
#else
#  include <deal.II/base/subscriptor.h>
#endif

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
MFOperator<dim, degree, number>::MFOperator(
  std::map<unsigned int, VariableAttributes>              _attributes_list,
  std::shared_ptr<const PDEOperator<dim, degree, number>> _pde_operator,
  Types::Index                                            _solve_block,
  Types::Index                                            _index,
  bool                                                    _use_local_mapping)
  : MATRIX_FREE_OPERATOR_BASE()
  , attributes_list(_attributes_list)
  , pde_operator(_pde_operator)
  , solve_block(_solve_block)
  , index(_index)
  , use_local_mapping(_use_local_mapping)
{}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::initialize(
  std::shared_ptr<const dealii::MatrixFree<dim, number, ScalarValue>> _data,
  const ElementVolume<dim, degree, number> &_element_volume_handler)
{
  data                   = _data;
  element_volume_handler = &_element_volume_handler;
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::types::global_dof_index
MFOperator<dim, degree, number>::m() const
{
  Assert(data.get() != nullptr, dealii::ExcNotInitialized());

  const unsigned int total_size =
    std::accumulate(selected_fields.begin(),
                    selected_fields.end(),
                    0U,
                    [this](unsigned int sum, unsigned int field)
                    {
                      return sum + data->get_vector_partitioner(field)->size();
                    });

  return total_size;
}

template <unsigned int dim, unsigned int degree, typename number>
number
MFOperator<dim, degree, number>::el([[maybe_unused]] const unsigned int &row,
                                    [[maybe_unused]] const unsigned int &col) const
{
  AssertThrow(false, FeatureNotImplemented("el()"));
  return 0.0;
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::clear()
{
  data.reset();
  diagonal_entries.reset();
  inverse_diagonal_entries.reset();
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::initialize_dof_vector(
  SolutionVector &dst,
  unsigned int    dof_handler_index) const
{
  data->initialize_dof_vector(dst, dof_handler_index);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::set_constrained_entries_to_one(SolutionVector &dst) const
{
  for (unsigned int j = 0; j < dealii::MFOperators::BlockHelper::n_blocks(dst); ++j)
    {
      const std::vector<unsigned int> &constrained_dofs =
        data->get_constrained_dofs(selected_fields[j]);
      for (const auto constrained_dof : constrained_dofs)
        {
          dealii::MFOperators::BlockHelper::subblock(dst, j).local_element(
            constrained_dof) = 1.0;
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
std::shared_ptr<const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>>
MFOperator<dim, degree, number>::get_matrix_free() const
{
  return data;
}

template <unsigned int dim, unsigned int degree, typename number>
const std::shared_ptr<
  dealii::DiagonalMatrix<typename MFOperator<dim, degree, number>::SolutionVector>> &
MFOperator<dim, degree, number>::get_matrix_diagonal_inverse() const
{
  Assert(inverse_diagonal_entries.get() != nullptr && inverse_diagonal_entries->m() > 0,
         dealii::ExcNotInitialized());
  return inverse_diagonal_entries;
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::compute_nonexplicit_auxiliary_update(
  std::vector<SolutionVector *>       &dst,
  const std::vector<SolutionVector *> &src) const
{
  Assert(!global_to_local_solution.empty(),
         dealii::ExcMessage(
           "The global to local solution mapping must not be empty. Make sure to call "
           "add_global_to_local_mapping() prior to any computations."));
  Assert(!dst.empty(), dealii::ExcMessage("The dst vector must not be empty"));
  Assert(!src.empty(), dealii::ExcMessage("The src vector must not be empty"));

  this->data->cell_loop(&MFOperator::compute_local_nonexplicit_auxiliary_update,
                        this,
                        dst,
                        src,
                        true);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::compute_residual(SolutionVector       &dst,
                                                  const SolutionVector &src) const
{
  Assert(!global_to_local_solution.empty(),
         dealii::ExcMessage(
           "The global to local solution mapping must not be empty. Make sure to call "
           "add_global_to_local_mapping() prior to any computations."));
  Assert(!src_solution_subset.empty(),
         dealii::ExcMessage("The src_solution_subset vector must not be empty"));
  Assert(dst.size() != 0,
         dealii::ExcMessage("The dst vector should not have size equal to 0"));
  Assert(src.size() != 0,
         dealii::ExcMessage("The src vector should not have size equal to 0"));

  this->data->cell_loop(&MFOperator::compute_local_residual, this, dst, src, true);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::vmult(SolutionVector       &dst,
                                       const SolutionVector &src) const
{
  Assert(!global_to_local_solution.empty(),
         dealii::ExcMessage(
           "The global to local solution mapping must not be empty. Make sure to call "
           "add_global_to_local_mapping() prior to any computations."));
  Assert(dst.size() != 0,
         dealii::ExcMessage("The dst vector should not have size equal to 0"));
  Assert(src.size() != 0,
         dealii::ExcMessage("The src vector should not have size equal to 0"));

  this->data->cell_loop(&MFOperator::compute_local_newton_update, this, dst, src, true);
}

// NOLINTBEGIN(readability-identifier-naming)

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::Tvmult(SolutionVector       &dst,
                                        const SolutionVector &src) const
{
  this->vmult(dst, src);
}

// NOLINTEND(readability-identifier-naming)

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::compute_local_rhs(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &_data,
  std::vector<SolutionVector *>                                          &dst,
  const std::vector<SolutionVector *>                                    &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  // The reason this is constructed here, rather than as a private member is because
  // compute_local_rhs is called by cell_loop, which multithreads. There would be data
  // races.
  FieldContainer<dim, degree, number> variable_list(_data,
                                                    attributes_list,
                                                    *element_volume_handler,
                                                    global_to_local_solution,
                                                    SolveType::ExplicitRHS);

  // Initialize, evaluate, and submit based on user function.
  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Grab the element volume
      const ScalarValue element_volume =
        variable_list.element_volume_handler->get_element_volume(cell);

      // Initialize, read DOFs, and set evaulation flags for each variable
      variable_list.reinit_and_eval(cell);

      // Evaluate at each quadrature point
      for (unsigned int quad = 0; quad < variable_list.get_n_q_points(); ++quad)
        {
          q_point = quad;
          pde_operator->compute_rhs(variable_list,
                                    variable_list.get_q_point_location(),
                                    element_volume);
        }

      // Integrate and add to global vector dst
      variable_list.integrate_and_distribute(dst);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::compute_local_newton_update(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &_data,
  SolutionVector                                                         &dst,
  const SolutionVector                                                   &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  FieldContainer<dim, degree, number> variable_list(_data,
                                                    attributes_list,
                                                    *element_volume_handler,
                                                    global_to_local_solution,
                                                    SolveType::NonexplicitLHS,
                                                    use_local_mapping);

  // Initialize, evaluate, and submit based on user function.
  variable_list.eval_local_operator(
    [this](FieldContainer<dim, degree, number>   &var_list,
           const dealii::Point<dim, ScalarValue> &q_point_loc,
           const ScalarValue                     &element_volume)
    {
      this->pde_operator->compute_nonexplicit_lhs(var_list,
                                                  q_point_loc,
                                                  element_volume,
                                                  solve_block,
                                                  index);
    },
    dst,
    src,
    src_solution_subset,
    cell_range);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::compute_diagonal(unsigned int field_index)
{
  inverse_diagonal_entries.reset(new dealii::DiagonalMatrix<SolutionVector>());
  SolutionVector &inverse_diagonal = inverse_diagonal_entries->get_vector();
  data->initialize_dof_vector(inverse_diagonal, field_index);
  const unsigned int dummy = 0;
  data->cell_loop(&MFOperator::local_compute_diagonal, this, inverse_diagonal, dummy);

  set_constrained_entries_to_one(inverse_diagonal);

  for (unsigned int i = 0; i < inverse_diagonal.locally_owned_size(); ++i)
    {
      Assert(inverse_diagonal.local_element(i) > static_cast<number>(0.0),
             dealii::ExcMessage(
               "No diagonal entry in a positive definite operator should be zero"));
      inverse_diagonal.local_element(i) =
        static_cast<number>(1.0) / inverse_diagonal.local_element(i);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
MFOperator<dim, degree, number>::local_compute_diagonal(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &_data,
  SolutionVector                                                         &dst,
  [[maybe_unused]] const unsigned int                                    &dummy,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  FieldContainer<dim, degree, number> variable_list(_data,
                                                    attributes_list,
                                                    *element_volume_handler,
                                                    global_to_local_solution,
                                                    SolveType::NonexplicitLHS,
                                                    use_local_mapping);

  for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
    {
      // Reinit the cell for all the dependencies
      variable_list.reinit_and_eval(cell);

      // To get the diagonal of the "matrix", repeatedly "multiply the matrix" by a test x
      // vector (change vector) and store the solution in the i-th position in the y
      // vector (the diagonal).
      // First (i=1) x vector: I 0 0 0
      // Second(i=2) x vector: 0 I 0 0 ...
      SolveGroup                        solve_group;
      GroupSolutionHandler<dim, number> group_solutions;

      // Submit zeros for everyting except the diagonals
      for (Types::FieldIndex field_index : solve_group.field_indices)
        {
          unsigned int n_dofs_per_cell = variable_list.get_dofs_per_component(
            field_index); // vector_feeval_ptr->dofs_per_component;
                          // scalar_feeval_ptr->dofs_per_cell;
          for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
            {
              int dof_index = i;
              for (unsigned int j = 0; j < n_dofs_per_cell; ++j)
                {
                  variable_list.set_dof_value(field_index,
                                              i == j ? Identity<Rank>() : Zero<Rank>(),
                                              j);
                }

              // Read plain dof values for non change src
              variable_list.read_dof_values(src_subset);

              // Evaluate the dependencies based on the flags
              variable_list.eval(global_variable_index);

              // Evaluate at each quadrature point
              for (unsigned int quad = 0; quad < get_n_q_points(); ++quad)
                {
                  q_point = quad;
                  compute_nonexplicit_lhs(*this, get_q_point_location(), element_volume);
                }

              // Integrate the diagonal
              integrate(global_variable_index);
              (*diagonal_ptr)[i] = feeval_ptr->get_dof_value(i);
            }
        }

      // Submit calculated diagonal values and distribute
      for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
        {
          if constexpr (std::is_same_v<DiagonalValueType, ScalarValue> || dim != 1)
            {
              feeval_ptr->submit_dof_value((*diagonal_ptr)[i], i);
            }
          else
            {
              feeval_ptr->submit_dof_value((*diagonal_ptr)[i][0], i);
            }
        }
      feeval_ptr->distribute_local_to_global(dst);
    }
}

#include "core/matrix_free_operator.inst"

PRISMS_PF_END_NAMESPACE
