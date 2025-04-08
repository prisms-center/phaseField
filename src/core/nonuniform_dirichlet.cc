#include <prismspf/core/nonuniform_dirichlet.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, fieldType field_type>
nonuniformDirichlet<dim, field_type>::nonuniformDirichlet(
  const unsigned int             &_index,
  const unsigned int             &_boundary_id,
  const userInputParameters<dim> &_user_inputs)
  : dealii::Function<dim>((field_type == fieldType::VECTOR) ? dim : 1)
  , index(_index)
  , boundary_id(_boundary_id)
  , user_inputs(&_user_inputs)
{}

// NOLINTBEGIN(readability-identifier-length)

template <int dim, fieldType field_type>
double
nonuniformDirichlet<dim, field_type>::value(
  const dealii::Point<dim>           &p,
  [[maybe_unused]] const unsigned int component) const
{
  // Initialize passed variables to zero
  double                 temp_scalar_value = 0.0;
  dealii::Vector<double> temp_vector_value(dim);

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

template <int dim, fieldType field_type>
void
nonuniformDirichlet<dim, field_type>::vector_value(const dealii::Point<dim> &p,
                                                   dealii::Vector<double>   &value) const
{
  // Initialize passed variables to zero
  double                 temp_scalar_value = 0.0;
  dealii::Vector<double> temp_vector_value(dim);

  // Pass variables to user-facing function to evaluate
  for (unsigned int i = 0; i < dim; i++)
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

template class nonuniformDirichlet<1, fieldType::SCALAR>;
template class nonuniformDirichlet<1, fieldType::VECTOR>;
template class nonuniformDirichlet<2, fieldType::SCALAR>;
template class nonuniformDirichlet<2, fieldType::VECTOR>;
template class nonuniformDirichlet<3, fieldType::SCALAR>;
template class nonuniformDirichlet<3, fieldType::VECTOR>;

PRISMS_PF_END_NAMESPACE