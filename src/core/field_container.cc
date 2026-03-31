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
  const DependencyMap                &dependency_map,
  const SolveBlock                   &_solve_group,
  const MatrixFree<dim, number>      &matrix_free)
  : field_attributes_ptr(&_field_attributes)
  , solution_indexer(&_solution_indexer)
  , solve_block(&_solve_group)
  , shared_feeval_scalar(matrix_free)
  , relative_level(_relative_level)
{
  const std::vector<FieldAttributes> &field_attributes = *field_attributes_ptr;

  // Initialize the feeval vectors
  feeval_deps_scalar.clear();
  feeval_deps_vector.clear();
  feeval_deps_scalar.resize(field_attributes.size());
  feeval_deps_vector.resize(field_attributes.size());

  // Loop over the dependency map
  for (const auto &[field_index, dependency] : dependency_map)
    {
      const auto mf_id_pair =
        solution_indexer->get_solution_level_and_block_index(field_index, relative_level);
      if (field_attributes[field_index].field_type == TensorRank::Scalar)
        {
          feeval_deps_scalar[field_index] = FEEValuationDeps<TensorRank::Scalar>(
            dependency,
            mf_id_pair,
            solve_block->id == solution_indexer->get_solve_group(field_index).id);
        }
      else if (field_attributes[field_index].field_type == TensorRank::Vector)
        {
          feeval_deps_vector[field_index] = FEEValuationDeps<TensorRank::Vector>(
            dependency,
            mf_id_pair,
            solve_block->id == solution_indexer->get_solve_group(field_index).id);
        }
      else
        {
          Assert(false, UnreachableCode());
        }
    }
}

#include "core/field_container.inst"

PRISMS_PF_END_NAMESPACE
