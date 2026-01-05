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

#include <functional>
#include <memory>
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
  const DependencySet                      &dependency_map,
  SolutionIndexer<dim, number>             &_solution_indexer,
  const ElementVolume<dim, degree, number> &_element_volume)
  : field_attributes_ptr(&_field_attributes)
  , solve_group(&_solve_group)
  , element_volume_handler(&_element_volume)
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
  /*else if vector */
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

// #include "core/field_container.inst"

PRISMS_PF_END_NAMESPACE
