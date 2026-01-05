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
  const dealii::MatrixFree<dim, number, dealii::VectorizedArray<number>> &data,
  const std::set<unsigned int>                                           &_field_indices,
  const std::vector<FieldAttributes> &_field_attributes,
  SolutionLevel<dim, number>         &_solution_level,
  const DependencySet                &dependency_map)
  : field_attributes_ptr(&_field_attributes)
{
  const std::vector<FieldAttributes> field_attributes;
  // Initialize the feeval vectors
  feeval_scalar.clear();
  feeval_vector.clear();
  feeval_scalar.resize(field_attributes.size());
  feeval_vector.resize(field_attributes.size());

  unsigned int block_index = 0;
  for (const auto &field_index : solve_group->field_indices)
    {
      const MatrixFree &matrix_free = solution_level.matrix_free;
      if (field_attributes[field_index].field_type == FieldInfo::TensorRank::Scalar)
        {
          feeval_scalar[field_index] = std::make_unique(matrix_free, block_index);
        }
      else if (field_attributes[field_index].field_type == FieldInfo::TensorRank::Vector)
        {
          feeval_vector[field_index] = std::make_unique(matrix_free, block_index);
        }
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
DSTContainer<dim, degree, number>::reinit(unsigned int cell)
{
  const DependencySet &dependency_map; // rhs or lhs
  for (const auto &[field_index, dependency] : dependency_map)
    {
      feeval_deps_scalar[field_index].reinit(cell);
    }
}

template <unsigned int dim, unsigned int degree, typename number>
void
DSTContainer<dim, degree, number>::integrate(unsigned int field_index,
                                             EvalFlags    submission_flags)
{
  const EvalFlags submission_flags = equation_type == EquationType::RHS
                                       ? field_attributes[field_index].eval_flags_rhs
                                       : field_attributes[field_index].eval_flags_lhs;
  if (field_attributes[field_index].field_type == TensorRank::Scalar)
    {
      feeval_scalar[field_index]->integrate(submission_flags);
    }
  else /* vector */
    {
      feeval_vector[field_index]->integrate(submission_flags);
    }

  /* AssertThrow(false,
              dealii::ExcMessage(
                "Integrate called for a solve type that is not NonexplicitLHS.")); */
}

template <unsigned int dim, unsigned int degree, typename number>
void
DSTContainer<dim, degree, number>::integrate_and_distribute(unsigned int field_index,
                                                            EvalFlags    submission_flags,
                                                            SolutionVector &dst)
{
  const std::vector<FieldAttributes> field_attributes = *field_attributes_ptr;

  if (field_attributes[field_index].field_type == TensorRank::Scalar)
    {
      feeval_scalar[field_index]->integrate_scatter(submission_flags, dst);
    }
  else /* vector */
    {
      feeval_vector[field_index]->integrate_scatter(submission_flags, dst);
    }
}

// #include "core/field_container.inst"

PRISMS_PF_END_NAMESPACE
