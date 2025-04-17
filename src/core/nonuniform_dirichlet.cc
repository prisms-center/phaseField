// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/nonuniform_dirichlet.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, typename number>
nonuniformDirichlet<dim, number>::nonuniformDirichlet(
  unsigned int                    _index,
  unsigned int                    _boundary_id,
  const userInputParameters<dim> &_user_inputs,
  unsigned int                    spacedim)
  : dealii::Function<dim, number>(spacedim)
  , index(_index)
  , boundary_id(_boundary_id)
  , user_inputs(&_user_inputs)
{}

// NOLINTBEGIN(readability-identifier-length)

template <int dim, typename number>
number
nonuniformDirichlet<dim, number>::value(
  const dealii::Point<dim>           &p,
  [[maybe_unused]] const unsigned int component) const
{
  // Initialize passed variables to zero
  number                 temp_scalar_value = 0.0;
  dealii::Vector<number> temp_vector_value(dim);

  // Pass variables to user-facing function to evaluate
  custom_nonuniform_dirichlet.set_nonuniform_dirichlet(index,
                                                       boundary_id,
                                                       0,
                                                       p,
                                                       temp_scalar_value,
                                                       temp_vector_value(0),
                                                       *user_inputs);

  return temp_scalar_value;
}

template <int dim, typename number>
void
nonuniformDirichlet<dim, number>::vector_value(const dealii::Point<dim> &p,
                                               dealii::Vector<number>   &value) const
{
  // Initialize passed variables to zero
  number                 temp_scalar_value = 0.0;
  dealii::Vector<number> temp_vector_value(dim);

  // Pass variables to user-facing function to evaluate
  for (int i = 0; i < dim; i++)
    {
      custom_nonuniform_dirichlet.set_nonuniform_dirichlet(index,
                                                           boundary_id,
                                                           i,
                                                           p,
                                                           temp_scalar_value,
                                                           temp_vector_value(i),
                                                           *user_inputs);
    }

  value = temp_vector_value;
}

// NOLINTEND(readability-identifier-length)

template class nonuniformDirichlet<1, double>;
template class nonuniformDirichlet<2, double>;
template class nonuniformDirichlet<3, double>;
template class nonuniformDirichlet<1, float>;
template class nonuniformDirichlet<2, float>;
template class nonuniformDirichlet<3, float>;

PRISMS_PF_END_NAMESPACE