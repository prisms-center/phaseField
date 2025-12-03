// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/types.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/exceptions.h>
#include <prismspf/core/matrix_free_operator.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>
#include <prismspf/core/variable_container.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/utilities/element_volume.h>

#include <prismspf/config.h>

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
MatrixFreeOperator<dim, degree, number>::MatrixFreeOperator(
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
MatrixFreeOperator<dim, degree, number>::initialize(
  std::shared_ptr<const dealii::MatrixFree<dim, number, SizeType>> _data,
  const ElementVolume<dim, degree, number> &_element_volume_handler,
  const std::vector<unsigned int>          &selected_field_indexes)
{
  data                   = _data;
  element_volume_handler = &_element_volume_handler;

  selected_fields.clear();
  if (selected_field_indexes.empty())
    {
      for (unsigned int i = 0; i < data->n_components(); ++i)
        {
          selected_fields.push_back(i);
        }
    }
  else
    {
      for (unsigned int i = 0; i < selected_field_indexes.size(); ++i)
        {
          AssertIndexRange(selected_field_indexes[i], data->n_components());
          for (unsigned int j = 0; j < selected_field_indexes.size(); ++j)
            {
              if (j != i)
                {
                  Assert(selected_field_indexes[j] != selected_field_indexes[i],
                         dealii::ExcMessage("Given row indices must be unique"));
                }
            }
          selected_fields.push_back(selected_field_indexes[i]);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
dealii::types::global_dof_index
MatrixFreeOperator<dim, degree, number>::m() const
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
MatrixFreeOperator<dim, degree, number>::el(
  [[maybe_unused]] const unsigned int &row,
  [[maybe_unused]] const unsigned int &col) const
{
  AssertThrow(false, FeatureNotImplemented("el()"));
  return 0.0;
}

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::clear()
{
  data.reset();
  diagonal_entries.reset();
  inverse_diagonal_entries.reset();
  global_to_local_solution.clear();
  src_solution_subset.clear();
  element_volume_handler = nullptr;
}

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::initialize_dof_vector(
  VectorType  &dst,
  unsigned int dof_handler_index) const
{
  data->initialize_dof_vector(dst, dof_handler_index);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::set_constrained_entries_to_one(
  VectorType &dst) const
{
  for (unsigned int j = 0; j < dealii::MatrixFreeOperators::BlockHelper::n_blocks(dst);
       ++j)
    {
      const std::vector<unsigned int> &constrained_dofs =
        data->get_constrained_dofs(selected_fields[j]);
      for (const auto constrained_dof : constrained_dofs)
        {
          dealii::MatrixFreeOperators::BlockHelper::subblock(dst, j).local_element(
            constrained_dof) = 1.0;
        }
    }
}

// cppcheck-suppress-begin passedByValue

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::add_global_to_local_mapping(
  std::vector<Types::Index> _global_to_local_solution)
{
  global_to_local_solution = _global_to_local_solution;
}

// cppcheck-suppress-end passedByValue

template <unsigned int dim, unsigned int degree, typename number>
std::shared_ptr<const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>>>
MatrixFreeOperator<dim, degree, number>::get_matrix_free() const
{
  return data;
}

template <unsigned int dim, unsigned int degree, typename number>
const std::shared_ptr<
  dealii::DiagonalMatrix<typename MatrixFreeOperator<dim, degree, number>::VectorType>> &
MatrixFreeOperator<dim, degree, number>::get_matrix_diagonal_inverse() const
{
  Assert(inverse_diagonal_entries.get() != nullptr && inverse_diagonal_entries->m() > 0,
         dealii::ExcNotInitialized());
  return inverse_diagonal_entries;
}

// cppcheck-suppress-begin passedByValue

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::add_src_solution_subset(
  std::vector<VectorType *> _src_solution_subset)
{
  src_solution_subset = _src_solution_subset;
}

// cppcheck-suppress-end passedByValue

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::compute_explicit_update(
  std::vector<VectorType *>       &dst,
  const std::vector<VectorType *> &src) const
{
  Assert(!global_to_local_solution.empty(),
         dealii::ExcMessage(
           "The global to local solution mapping must not be empty. Make sure to call "
           "add_global_to_local_mapping() prior to any computations."));
  Assert(!dst.empty(), dealii::ExcMessage("The dst vector must not be empty"));
  Assert(!src.empty(), dealii::ExcMessage("The src vector must not be empty"));

  this->data->cell_loop(&MatrixFreeOperator::compute_local_explicit_update,
                        this,
                        dst,
                        src,
                        true);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::compute_postprocess_explicit_update(
  std::vector<VectorType *>       &dst,
  const std::vector<VectorType *> &src) const
{
  Assert(!global_to_local_solution.empty(),
         dealii::ExcMessage(
           "The global to local solution mapping must not be empty. Make sure to call "
           "add_global_to_local_mapping() prior to any computations."));
  Assert(!dst.empty(), dealii::ExcMessage("The dst vector must not be empty"));
  Assert(!src.empty(), dealii::ExcMessage("The src vector must not be empty"));

  this->data->cell_loop(&MatrixFreeOperator::compute_local_postprocess_explicit_update,
                        this,
                        dst,
                        src,
                        true);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::compute_nonexplicit_auxiliary_update(
  std::vector<VectorType *>       &dst,
  const std::vector<VectorType *> &src) const
{
  Assert(!global_to_local_solution.empty(),
         dealii::ExcMessage(
           "The global to local solution mapping must not be empty. Make sure to call "
           "add_global_to_local_mapping() prior to any computations."));
  Assert(!dst.empty(), dealii::ExcMessage("The dst vector must not be empty"));
  Assert(!src.empty(), dealii::ExcMessage("The src vector must not be empty"));

  this->data->cell_loop(&MatrixFreeOperator::compute_local_nonexplicit_auxiliary_update,
                        this,
                        dst,
                        src,
                        true);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::compute_residual(VectorType       &dst,
                                                          const VectorType &src) const
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

  this->data->cell_loop(&MatrixFreeOperator::compute_local_residual,
                        this,
                        dst,
                        src,
                        true);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::vmult(VectorType       &dst,
                                               const VectorType &src) const
{
  Assert(!global_to_local_solution.empty(),
         dealii::ExcMessage(
           "The global to local solution mapping must not be empty. Make sure to call "
           "add_global_to_local_mapping() prior to any computations."));
  Assert(dst.size() != 0,
         dealii::ExcMessage("The dst vector should not have size equal to 0"));
  Assert(src.size() != 0,
         dealii::ExcMessage("The src vector should not have size equal to 0"));

  this->data->cell_loop(&MatrixFreeOperator::compute_local_newton_update,
                        this,
                        dst,
                        src,
                        true);
}

// NOLINTBEGIN(readability-identifier-naming)

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::Tvmult(VectorType       &dst,
                                                const VectorType &src) const
{
  this->vmult(dst, src);
}

// NOLINTEND(readability-identifier-naming)

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::compute_local_explicit_update(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &_data,
  std::vector<VectorType *>                                              &dst,
  const std::vector<VectorType *>                                        &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  VariableContainer<dim, degree, number> variable_list(_data,
                                                       attributes_list,
                                                       *element_volume_handler,
                                                       global_to_local_solution,
                                                       SolveType::ExplicitRHS);

  // Initialize, evaluate, and submit based on user function.
  variable_list.eval_local_operator(
    [this](VariableContainer<dim, degree, number> &var_list,
           const dealii::Point<dim, SizeType>     &q_point_loc,
           const SizeType                         &element_volume)
    {
      this->pde_operator->compute_explicit_rhs(var_list,
                                               q_point_loc,
                                               element_volume,
                                               solve_block);
    },
    dst,
    src,
    cell_range);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::compute_local_postprocess_explicit_update(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &_data,
  std::vector<VectorType *>                                              &dst,
  const std::vector<VectorType *>                                        &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  VariableContainer<dim, degree, number> variable_list(_data,
                                                       attributes_list,
                                                       *element_volume_handler,
                                                       global_to_local_solution,
                                                       SolveType::Postprocess);

  // Initialize, evaluate, and submit based on user function.
  variable_list.eval_local_operator(
    [this](VariableContainer<dim, degree, number> &var_list,
           const dealii::Point<dim, SizeType>     &q_point_loc,
           const SizeType                         &element_volume)
    {
      this->pde_operator->compute_postprocess_explicit_rhs(var_list,
                                                           q_point_loc,
                                                           element_volume,
                                                           solve_block);
    },
    dst,
    src,
    cell_range);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::compute_local_nonexplicit_auxiliary_update(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &_data,
  std::vector<VectorType *>                                              &dst,
  const std::vector<VectorType *>                                        &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  VariableContainer<dim, degree, number> variable_list(_data,
                                                       attributes_list,
                                                       *element_volume_handler,
                                                       global_to_local_solution,
                                                       SolveType::NonexplicitRHS);

  // Initialize, evaluate, and submit based on user function.
  variable_list.eval_local_operator(
    [this](VariableContainer<dim, degree, number> &var_list,
           const dealii::Point<dim, SizeType>     &q_point_loc,
           const SizeType                         &element_volume)
    {
      this->pde_operator->compute_nonexplicit_rhs(var_list,
                                                  q_point_loc,
                                                  element_volume,
                                                  solve_block,
                                                  index);
    },
    dst,
    src,
    cell_range);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::compute_local_residual(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &_data,
  VectorType                                                             &dst,
  [[maybe_unused]] const VectorType                                      &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  VariableContainer<dim, degree, number> variable_list(_data,
                                                       attributes_list,
                                                       *element_volume_handler,
                                                       global_to_local_solution,
                                                       SolveType::NonexplicitRHS);

  // Initialize, evaluate, and submit based on user function.
  variable_list.eval_local_operator(
    [this](VariableContainer<dim, degree, number> &var_list,
           const dealii::Point<dim, SizeType>     &q_point_loc,
           const SizeType                         &element_volume)
    {
      this->pde_operator->compute_nonexplicit_rhs(var_list,
                                                  q_point_loc,
                                                  element_volume,
                                                  solve_block,
                                                  index);
    },
    dst,
    src_solution_subset,
    cell_range);
}

template <unsigned int dim, unsigned int degree, typename number>
void
MatrixFreeOperator<dim, degree, number>::compute_local_newton_update(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &_data,
  VectorType                                                             &dst,
  const VectorType                                                       &src,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  VariableContainer<dim, degree, number> variable_list(_data,
                                                       attributes_list,
                                                       *element_volume_handler,
                                                       global_to_local_solution,
                                                       SolveType::NonexplicitLHS,
                                                       use_local_mapping);

  // Initialize, evaluate, and submit based on user function. Note that the src solution
  // subset must not include the src vector.
  variable_list.eval_local_operator(
    [this](VariableContainer<dim, degree, number> &var_list,
           const dealii::Point<dim, SizeType>     &q_point_loc,
           const SizeType                         &element_volume)
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
MatrixFreeOperator<dim, degree, number>::compute_diagonal(unsigned int field_index)
{
  inverse_diagonal_entries.reset(new dealii::DiagonalMatrix<VectorType>());
  VectorType &inverse_diagonal = inverse_diagonal_entries->get_vector();
  data->initialize_dof_vector(inverse_diagonal, field_index);
  const unsigned int dummy = 0;
  data->cell_loop(&MatrixFreeOperator::local_compute_diagonal,
                  this,
                  inverse_diagonal,
                  dummy);

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
MatrixFreeOperator<dim, degree, number>::local_compute_diagonal(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &_data,
  VectorType                                                             &dst,
  [[maybe_unused]] const unsigned int                                    &dummy,
  const std::pair<unsigned int, unsigned int> &cell_range) const
{
  // Constructor for FEEvaluation objects
  VariableContainer<dim, degree, number> variable_list(_data,
                                                       attributes_list,
                                                       *element_volume_handler,
                                                       global_to_local_solution,
                                                       SolveType::NonexplicitLHS,
                                                       use_local_mapping);

  // Initialize, evaluate, and submit diagonal based on user function.
  variable_list.eval_local_diagonal(
    [this](VariableContainer<dim, degree, number> &var_list,
           const dealii::Point<dim, SizeType>     &q_point_loc,
           const SizeType                         &element_volume)
    {
      this->pde_operator->compute_nonexplicit_lhs(var_list,
                                                  q_point_loc,
                                                  element_volume,
                                                  solve_block,
                                                  index);
    },
    dst,
    src_solution_subset,
    cell_range);
}

#include "core/matrix_free_operator.inst"

PRISMS_PF_END_NAMESPACE
