// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/pde_operator.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/config.h>
#include <prismspf/field_input/read_field_factory.h>

#include <memory>
#include <string>
#include <utility>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
InitialCondition<dim, degree, number>::InitialCondition(
  const unsigned int                                            &_index,
  const FieldInfo::TensorRank                                   &_field_type,
  const std::shared_ptr<const PDEOperator<dim, degree, number>> &_pde_operator)
  : dealii::Function<dim, number>((_field_type == FieldInfo::TensorRank::Vector) ? dim
                                                                                 : 1)
  , index(_index)
  , field_type(_field_type)
  , pde_operator(_pde_operator)
{}

// NOLINTBEGIN(readability-identifier-length)

template <unsigned int dim, unsigned int degree, typename number>
void
InitialCondition<dim, degree, number>::vector_value(const dealii::Point<dim> &p,
                                                    dealii::Vector<number>   &value) const
{
  // Initialize passed variables to zero
  dealii::Vector<number> vector_value(dim);

  // Pass variables to user-facing function to evaluate
  // TODO (landinjm): For the values that don't make sense we should pass invalid values
  // (e.g., underflow values) so identity the subtle bug. For example, the vector values
  // and component in this scalar condition should be uint underflow and NaN.
  if (field_type == FieldInfo::TensorRank::Scalar)
    {
      pde_operator->set_initial_condition(index, 0, p, vector_value[0], vector_value[0]);
    }
  else if (field_type == FieldInfo::TensorRank::Vector)
    {
      for (unsigned int i = 0; i < dim; i++)
        {
          pde_operator->set_initial_condition(index,
                                              i,
                                              p,
                                              vector_value[0],
                                              vector_value[i]);
        }
    }

  value = vector_value;
}

// NOLINTEND(readability-identifier-length)

template <unsigned int dim, typename number>
ReadInitialCondition<dim, number>::ReadInitialCondition(
  std::string                       _field_name,
  const FieldInfo::TensorRank      &_field_type,
  const InitialConditionFile       &ic_file,
  const SpatialDiscretization<dim> &spatial_discretization)
  : dealii::Function<dim, number>((_field_type == FieldInfo::TensorRank::Vector) ? dim
                                                                                 : 1)
  , field_name(_field_name)
  , field_type(_field_type)
  , reader(create_reader<dim, number>(ic_file, spatial_discretization))
{}

// NOLINTBEGIN(readability-identifier-length)

template <unsigned int dim, typename number>
void
ReadInitialCondition<dim, number>::vector_value(const dealii::Point<dim> &p,
                                                dealii::Vector<number>   &value) const
{
  // Initialize passed variables to zero
  dealii::Vector<number> vector_value(dim);

  // Get the value of the field at the point
  if (field_type == FieldInfo::TensorRank::Scalar)
    {
      vector_value[0] = reader->get_scalar_value(p, field_name);
    }
  else if (field_type == FieldInfo::TensorRank::Vector)
    {
      vector_value = reader->get_vector_value(p, field_name);
    }

  value = vector_value;
}

// NOLINTEND(readability-identifier-length)

#include "core/initial_conditions.inst"

PRISMS_PF_END_NAMESPACE
