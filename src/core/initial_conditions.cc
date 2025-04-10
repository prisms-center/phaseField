#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/core/initial_conditions.h>
#include <prismspf/core/type_enums.h>

#include <prismspf/user_inputs/user_input_parameters.h>

#include <prismspf/config.h>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim>
initialCondition<dim>::initialCondition(const unsigned int             &_index,
                                        const fieldType                &field_type,
                                        const userInputParameters<dim> &_user_inputs)
  : dealii::Function<dim>((field_type == fieldType::VECTOR) ? dim : 1)
  , index(_index)
  , user_inputs(&_user_inputs)
{}

// NOLINTBEGIN(readability-identifier-length)

template <int dim>
void
initialCondition<dim>::vector_value(const dealii::Point<dim> &p,
                                    dealii::Vector<double>   &value) const
{
  // Initialize passed variables to zero
  dealii::Vector<double> vector_value(dim);

  // Pass variables to user-facing function to evaluate
  for (unsigned int i = 0; i < dim; i++)
    {
      custom_initial_condition.set_initial_condition(index,
                                                     i,
                                                     p,
                                                     vector_value(0),
                                                     vector_value(i),
                                                     *user_inputs);
    }

  value = vector_value;
}

// NOLINTEND(readability-identifier-length)

INSTANTIATE_UNI_TEMPLATE(initialCondition)

PRISMS_PF_END_NAMESPACE