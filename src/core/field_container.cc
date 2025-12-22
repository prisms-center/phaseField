// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/exceptions.h>
#include <prismspf/core/field_container.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/utilities/element_volume.h>

#include <prismspf/config.h>

#include "prismspf/core/dependencies.h"
#include "prismspf/core/field_attributes.h"
#include "prismspf/core/solution_indexer.h"
#include "prismspf/core/solve_group.h"

#include <algorithm>
#include <functional>
#include <map>
#include <memory>
#include <ranges>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
FieldContainer<dim, degree, number>::FieldContainer(
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
  const std::vector<FieldAttributes>       &_field_attributes,
  const SolveGroup                         &_solve_group,
  SolutionIndexer<dim, number>             &_solution_indexer,
  const ElementVolume<dim, degree, number> &_element_volume,
  const EquationType                        _equation_type,
  bool                                      use_local_mapping)
  : field_attributes_ptr(&_field_attributes)
  , solve_group(&_solve_group)
  , element_volume_handler(&_element_volume)
  , equation_type(equation_type)
{
  const std::vector<FieldAttributes> field_attributes;
  // Initialize the feeval vectors
  feeval_deps_scalar.clear();
  feeval_deps_vector.clear();
  dst_feeval_scalar.clear();
  dst_feeval_vector.clear();
  feeval_deps_scalar.resize(field_attributes.size());
  feeval_deps_vector.resize(field_attributes.size());
  dst_feeval_scalar.resize(field_attributes.size());
  dst_feeval_vector.resize(field_attributes.size());

  const DependencySet &dependency_map = equation_type == EquationType::RHS
                                          ? solve_group->dependencies_rhs
                                          : solve_group->dependencies_lhs;
  for (const auto &[field_index, dependency] : dependency_map)
    {
      const auto mf_id_pair =
        solution_indexer->get_solution_level_and_block_index(field_index, relative_level);
      if (field_attributes[field_index].field_type == FieldInfo::TensorRank::Scalar)
        {
          feeval_deps_scalar[field_index] = {dependency, mf_id_pair};
        }
      else if (field_attributes[field_index].field_type == FieldInfo::TensorRank::Vector)
        {
          feeval_deps_vector[field_index] = {dependency, mf_id_pair};
        }
      // else Assert unreachable
    }
  MatrixFree   matrix_free;
  unsigned int block_index = 0;
  for (const auto &field_index : solve_group->field_indices)
    {
      if (field_attributes[field_index].field_type == FieldInfo::TensorRank::Scalar)
        {
          dst_feeval_scalar[field_index] = std::make_unique(matrix_free, block_index);
        }
      else if (field_attributes[field_index].field_type == FieldInfo::TensorRank::Vector)
        {
          dst_feeval_vector[field_index] = std::make_unique(matrix_free, block_index);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::eval_local_diagonal(
  const std::function<void(FieldContainer &,
                           const dealii::Point<dim, ScalarValue> &,
                           const ScalarValue &)> &func,
  SolutionVector                                 &dst,
  const std::vector<SolutionVector *>            &src_subset,
  const std::pair<unsigned int, unsigned int>    &cell_range)
{
  Assert(subset_attributes->size() == 1,
         dealii::ExcMessage(
           "For nonexplicit solves, subset attributes should only be 1 variable."));

  const auto &global_var_index = subset_attributes->begin()->first;
  const auto &field_type = subset_attributes->begin()->second.field_info.tensor_rank;
  feevaluation_exists(global_var_index, DependencyType::Change);
  auto &feeval_variant = feeval_vector[(global_var_index * max_dependency_types) +
                                       static_cast<Types::Index>(DependencyType::Change)];

  auto process_feeval = [&](auto &feeval_ptr, auto &diag_ptr)
  {
    using FEEvaluationType = std::decay_t<decltype(*feeval_ptr)>;
    using DiagonalType     = std::decay_t<decltype(*diag_ptr)>;

    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        eval_cell_diagonal<FEEvaluationType, DiagonalType>(feeval_ptr,
                                                           diag_ptr,
                                                           cell,
                                                           global_var_index,
                                                           func,
                                                           dst,
                                                           src_subset);
      }
  };

  if (field_type == FieldInfo::TensorRank::Scalar)
    {
      ScalarFEEvaluation *scalar_feeval_ptr = nullptr;

      if constexpr (dim == 1)
        {
          scalar_feeval_ptr = feeval_variant.get();
        }
      else
        {
          scalar_feeval_ptr = extract_feeval_ptr<ScalarFEEvaluation>(feeval_variant);
        }

      n_dofs_per_cell       = scalar_feeval_ptr->dofs_per_cell;
      diagonal              = std::make_unique<ScalarDiagonal>(n_dofs_per_cell);
      auto *scalar_diag_ptr = extract_diagonal_ptr<ScalarDiagonal>(diagonal);
      process_feeval(scalar_feeval_ptr, scalar_diag_ptr);
    }
  else if (field_type == FieldInfo::TensorRank::Vector)
    {
      VectorFEEvaluation *vector_feeval_ptr = nullptr;

      if constexpr (dim == 1)
        {
          vector_feeval_ptr = feeval_variant.get();
        }
      else
        {
          vector_feeval_ptr = extract_feeval_ptr<VectorFEEvaluation>(feeval_variant);
        }

      n_dofs_per_cell       = vector_feeval_ptr->dofs_per_component;
      diagonal              = std::make_unique<VectorDiagonal>(n_dofs_per_cell);
      auto *vector_diag_ptr = extract_diagonal_ptr<VectorDiagonal>(diagonal);
      process_feeval(vector_feeval_ptr, vector_diag_ptr);
    }
  else
    {
      Assert(false, UnreachableCode());
    }
}

template <unsigned int dim, unsigned int degree, typename number>
unsigned int
FieldContainer<dim, degree, number>::get_n_q_points() const
{}

template <unsigned int dim, unsigned int degree, typename number>
dealii::Point<dim, typename FieldContainer<dim, degree, number>::ScalarValue>
FieldContainer<dim, degree, number>::get_q_point_location() const
{
  Types::Index                        field_index = *(solve_group->field_indices.begin());
  const std::vector<FieldAttributes> &field_attributes = *field_attributes_ptr;
  if (field_attributes[field_index].field_type == TensorRank::Scalar)
    {
      return dst_feeval_scalar[field_index]->quadrature_point(q_point);
    }
  else /* vector */
    {
      return dst_feeval_vector[field_index]->quadrature_point(q_point);
    }
  /* unreachable */
  return dealii::Point<dim, ScalarValue>();
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::reinit(unsigned int cell)
{
  const DependencySet &dependency_map; // rhs or lhs
  for (const auto &[field_index, dependency] : dependency_map)
    {
      feeval_deps_scalar[field_index].reinit(cell);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::eval()
{
  const DependencySet &dependency_map; // rhs or lhs
  for (const auto &[field_index, dependency] : dependency_map)
    {
      feeval_deps_scalar[field_index].eval();
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::reinit_and_eval(unsigned int cell)
{
  const DependencySet &dependency_map; // rhs or lhs
  for (const auto &[field_index, dependency] : dependency_map)
    {
      feeval_deps_scalar[field_index].reinit_and_eval(cell);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::integrate()
{
  const std::vector<FieldAttributes> field_attributes = *field_attributes_ptr;
  for (const Types::Index &field_index : solve_group->field_indices)
    {
      const EvalFlags submission_flags = equation_type == EquationType::RHS
                                           ? field_attributes[field_index].eval_flags_rhs
                                           : field_attributes[field_index].eval_flags_lhs;
      if (field_attributes[field_index].field_type == TensorRank::Scalar)
        {
          dst_feeval_scalar[field_index]->integrate(submission_flags);
        }
      else /* vector */
        {
          dst_feeval_vector[field_index]->integrate(submission_flags);
        }
    }
  /* AssertThrow(false,
              dealii::ExcMessage(
                "Integrate called for a solve type that is not NonexplicitLHS.")); */
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::integrate_and_distribute()
{
  const std::vector<FieldAttributes> field_attributes = *field_attributes_ptr;
  for (const Types::Index &field_index : solve_group->field_indices)
    {
      const EvalFlags submission_flags = equation_type == EquationType::RHS
                                           ? field_attributes[field_index].eval_flags_rhs
                                           : field_attributes[field_index].eval_flags_lhs;
      SolutionVector &dst              = equation_type == EquationType::RHS
                                           ? solutions.new_solutions(field_index, relative_level)
                                           : solutions.change_solutions(field_index, relative_level);
      ;
      if (field_attributes[field_index].field_type == TensorRank::Scalar)
        {
          dst_feeval_scalar[field_index]->integrate_scatter(submission_flags, dst);
        }
      else /* vector */
        {
          dst_feeval_vector[field_index]->integrate_scatter(submission_flags, dst);
        }
    }
}

#include "core/variable_container.inst"

PRISMS_PF_END_NAMESPACE
