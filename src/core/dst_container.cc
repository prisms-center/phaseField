// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/vectorization.h>
#include <deal.II/matrix_free/evaluation_flags.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/dependencies.h>
#include <prismspf/core/dst_container.h>
#include <prismspf/core/exceptions.h>
#include <prismspf/core/field_attributes.h>
#include <prismspf/core/solution_indexer.h>
#include <prismspf/core/solve_group.h>
#include <prismspf/core/type_enums.h>
#include <prismspf/core/types.h>
#include <prismspf/core/variable_attributes.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
DSTContainer<dim, degree, number>::DSTContainer(
  const std::set<unsigned int>       &_field_indices,
  const std::vector<FieldAttributes> &_field_attributes,
  const MatrixFree                   &matrix_free,
  const std::vector<unsigned int>    &_field_to_block_index)
  : field_attributes_ptr(&_field_attributes)
  , field_indices(&_field_indices)
{
  const std::vector<FieldAttributes> field_attributes;
  // Initialize the feeval vectors
  feeval_scalar.clear();
  feeval_vector.clear();
  feeval_scalar.resize(field_attributes.size());
  feeval_vector.resize(field_attributes.size());
  integration_flags = std::vector<EvalFlags>(field_attributes.size(), EvalFlags::nothing);
  for (const auto &field_index : *field_indices)
    {
      const unsigned int block_index = _field_to_block_index[field_index];
      if (field_attributes[field_index].field_type == FieldInfo::TensorRank::Scalar)
        {
          feeval_scalar[field_index] = std::make_unique(matrix_free, block_index);
        }
      else
        {
          feeval_vector[field_index] = std::make_unique(matrix_free, block_index);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
DSTContainer<dim, degree, number>::reinit(unsigned int cell)
{
  for (const unsigned int &field_index : *field_indices)
    {
      if ((*field_attributes_ptr)[field_index].field_type == TensorRank::Scalar)
        {
          feeval_scalar[field_index].reinit(cell);
        }
      else
        {
          feeval_vector[field_index].reinit(cell);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
DSTContainer<dim, degree, number>::integrate(unsigned int field_index)
{
  if ((*field_attributes_ptr)[field_index].field_type == TensorRank::Scalar)
    {
      feeval_scalar[field_index]->integrate(integration_flags[field_index]);
    }
  else /* vector */
    {
      feeval_vector[field_index]->integrate(integration_flags[field_index]);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
DSTContainer<dim, degree, number>::integrate_and_distribute(unsigned int    field_index,
                                                            SolutionVector &dst)
{
  const std::vector<FieldAttributes> field_attributes = *field_attributes_ptr;

  if (field_attributes[field_index].field_type == TensorRank::Scalar)
    {
      feeval_scalar[field_index]->integrate_scatter(integration_flags[field_index], dst);
    }
  else /* vector */
    {
      feeval_vector[field_index]->integrate_scatter(integration_flags[field_index], dst);
    }
}

// #include "core/field_container.inst"

PRISMS_PF_END_NAMESPACE
