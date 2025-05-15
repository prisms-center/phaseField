// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/nonuniform_dirichlet.h>
#include <prismspf/core/pde_operator.h>

#include <prismspf/config.h>

#include <memory>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
nonuniformDirichlet<dim, degree, number>::nonuniformDirichlet(
  unsigned int                                                   _index,
  unsigned int                                                   _boundary_id,
  const std::shared_ptr<const PDEOperator<dim, degree, number>> &_pde_operator,
  unsigned int                                                   spacedim)
  : dealii::Function<dim, number>(spacedim)
  , index(_index)
  , boundary_id(_boundary_id)
  , pde_operator(_pde_operator)
{}

// NOLINTBEGIN(readability-identifier-length)

template <unsigned int dim, unsigned int degree, typename number>
number
nonuniformDirichlet<dim, degree, number>::value(
  const dealii::Point<dim>           &p,
  [[maybe_unused]] const unsigned int component) const
{
  // Initialize passed variables to zero
  number                 temp_scalar_value = 0.0;
  dealii::Vector<number> temp_vector_value(dim);

  // Pass variables to user-facing function to evaluate
  pde_operator->set_nonuniform_dirichlet(index,
                                         boundary_id,
                                         0,
                                         p,
                                         temp_scalar_value,
                                         temp_vector_value[0]);

  return temp_scalar_value;
}

template <unsigned int dim, unsigned int degree, typename number>
void
nonuniformDirichlet<dim, degree, number>::vector_value(
  const dealii::Point<dim> &p,
  dealii::Vector<number>   &value) const
{
  // TODO (landinjm): I think this function is not called for 1D vector and might break
  // when the user goes from 2D to 1D vector fields.

  // Initialize passed variables to zero
  number                 temp_scalar_value = 0.0;
  dealii::Vector<number> temp_vector_value(dim);

  // Pass variables to user-facing function to evaluate
  for (unsigned int i = 0; i < dim; i++)
    {
      pde_operator->set_nonuniform_dirichlet(index,
                                             boundary_id,
                                             0,
                                             p,
                                             temp_scalar_value,
                                             temp_vector_value[i]);
    }

  value = temp_vector_value;
}

// NOLINTEND(readability-identifier-length)

INSTANTIATE_TRI_TEMPLATE(nonuniformDirichlet)

PRISMS_PF_END_NAMESPACE