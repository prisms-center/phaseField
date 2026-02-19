// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>

#include <prismspf/core/exceptions.h>
#include <prismspf/core/field_container.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
FieldContainer<dim, degree, number>::FieldContainer(
  const std::vector<FieldAttributes> &_field_attributes,
  const SolutionIndexer<dim, number> &_solution_indexer,
  unsigned int                        _relative_level,
  const DependencySet                &dependency_map,
  const SolveGroup                   &_solve_group,
  const MatrixFree                   &matrix_free)
  : field_attributes_ptr(&_field_attributes)
  , solution_indexer(&_solution_indexer)
  , relative_level(_relative_level)
  , solve_group(&_solve_group)
  , shared_feeval_scalar(matrix_free)
{
  const std::vector<FieldAttributes> &field_attributes = *field_attributes_ptr;
  // Initialize the feeval vectors
  feeval_deps_scalar.clear();
  feeval_deps_vector.clear();
  feeval_deps_scalar.resize(field_attributes.size());
  feeval_deps_vector.resize(field_attributes.size());

  for (const auto &[field_index, dependency] : dependency_map)
    {
      const auto mf_id_pair =
        solution_indexer->get_solution_level_and_block_index(field_index, relative_level);
      if (field_attributes[field_index].field_type == FieldInfo::TensorRank::Scalar)
        {
          feeval_deps_scalar[field_index] = FEEValuationDeps<ScalarFEEvaluation>(
            dependency,
            mf_id_pair,
            solve_group->id == solution_indexer->get_solve_group(field_index).id);
        }
      else if (field_attributes[field_index].field_type == FieldInfo::TensorRank::Vector)
        {
          feeval_deps_vector[field_index] = FEEValuationDeps<VectorFEEvaluation>(
            dependency,
            mf_id_pair,
            solve_group->id == solution_indexer->get_solve_group(field_index).id);
        }
      // else Assert unreachable
    }
  //================================================================================
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::reinit(unsigned int cell)
{
  for (auto &fe_eval : feeval_deps_scalar)
    {
      fe_eval.reinit(cell);
    }
  for (auto &fe_eval : feeval_deps_vector)
    {
      fe_eval.reinit(cell);
    }
  //================================================================================
  shared_feeval_scalar.reinit(cell);
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::eval(const BlockVector *src_solutions)
{
  for (auto &fe_eval : feeval_deps_scalar)
    {
      fe_eval.eval(src_solutions);
    }
  for (auto &fe_eval : feeval_deps_vector)
    {
      fe_eval.eval(src_solutions);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::reinit_and_eval(unsigned int       cell,
                                                     const BlockVector *src_solutions)
{
  for (auto &fe_eval : feeval_deps_scalar)
    {
      fe_eval.reinit_and_eval(cell, src_solutions);
    }
  for (auto &fe_eval : feeval_deps_vector)
    {
      fe_eval.reinit_and_eval(cell, src_solutions);
    }
  //================================================================================
  shared_feeval_scalar.reinit(cell);
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::integrate()
{
  const std::vector<FieldAttributes> &field_attributes = *field_attributes_ptr;
  for (const Types::Index &field_index : solve_group->field_indices)
    {
      if (field_attributes[field_index].field_type == TensorRank::Scalar)
        {
          feeval_deps_scalar[field_index].integrate();
        }
      else /* vector */
        {
          feeval_deps_vector[field_index].integrate();
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::distribute(BlockVector *dst_solutions)
{
  const std::vector<FieldAttributes> &field_attributes = *field_attributes_ptr;
  for (const Types::Index &field_index : solve_group->field_indices)
    {
      if (field_attributes[field_index].field_type == TensorRank::Scalar)
        {
          feeval_deps_scalar[field_index].distribute(dst_solutions);
        }
      else /* vector */
        {
          feeval_deps_vector[field_index].distribute(dst_solutions);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::integrate_and_distribute(BlockVector *dst_solutions)
{
  const std::vector<FieldAttributes> &field_attributes = *field_attributes_ptr;
  for (const Types::Index &field_index : solve_group->field_indices)
    {
      if (field_attributes[field_index].field_type == TensorRank::Scalar)
        {
          feeval_deps_scalar[field_index].integrate_and_distribute(dst_solutions);
        }
      else /* vector */
        {
          feeval_deps_vector[field_index].integrate_and_distribute(dst_solutions);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::feevaluation_size_valid(
  Types::Index field_index) const
{
  // TODO
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::feevaluation_exists(
  Types::Index field_index,
  Types::Index dependency_index) const
{
  // TODO
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::access_valid(
  Types::Index                             field_index,
  DependencyType                           dependency_type,
  dealii::EvaluationFlags::EvaluationFlags flag) const
{
  // TODO
}

template <unsigned int dim, unsigned int degree, typename number>
void
FieldContainer<dim, degree, number>::submission_valid(
  Types::Index   field_index,
  DependencyType dependency_type) const
{
  // TODO
}

#include "core/field_container.inst"

PRISMS_PF_END_NAMESPACE
