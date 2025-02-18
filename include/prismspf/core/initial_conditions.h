#ifndef initial_conditions_h
#define initial_conditions_h

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>
#include <deal.II/lac/vector.h>

#include <prismspf/config.h>
#include <prismspf/core/type_enums.h>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * \brief Forward declaration of user-facing implementation
 */
template <int dim>
class customInitialCondition;

/**
 * \brief Function for user-implemented initial conditions. These are only ever calculated
 * for explicit time dependent fields and implicit time dependent, as all others are
 * calculated at runtime.
 */
template <int dim>
class initialCondition : public dealii::Function<dim, double>
{
public:
  /**
   * \brief Constructor.
   */
  initialCondition(const unsigned int &_index, const fieldType &field_type);

  /**
   * \brief Scalar/Vector value.
   */
  void
  vector_value(const dealii::Point<dim> &p, dealii::Vector<double> &value) const override;

private:
  const unsigned int index;

  customInitialCondition<dim> custom_initial_condition;
};

template <int dim>
initialCondition<dim>::initialCondition(const unsigned int &_index,
                                        const fieldType    &field_type)
  : dealii::Function<dim>((field_type == fieldType::VECTOR) ? dim : 1)
  , index(_index)
{}

template <int dim>
inline void
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
                                                     vector_value(i));
    }

  value = vector_value;
}

/**
 * \brief User-facing implementation of initial conditions
 */
template <int dim>
class customInitialCondition
{
public:
  /**
   * \brief Constructor.
   */
  customInitialCondition() = default;

  /**
   * \brief Function that passes the value/vector and point that are set in the initial
   * condition.
   */
  void
  set_initial_condition(const unsigned int       &index,
                        const unsigned int       &component,
                        const dealii::Point<dim> &point,
                        double                   &scalar_value,
                        double                   &vector_component_value) const;
};

PRISMS_PF_END_NAMESPACE

#endif