#include <deal.II/base/exceptions.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <prismspf/core/type_enums.h>

#include <prismspf/utilities/compute_integral.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <int dim, int degree, typename number>
void
computeIntegral<dim, degree, number>::compute_integral(
  number                        &integral_value,
  const dealii::DoFHandler<dim> &dof_handler,
  const VectorType              &vector) const
{
  [[maybe_unused]] const unsigned int expected_components = 1;
  Assert(dof_handler.get_fe().n_components() == expected_components,
         dealii::ExcMessage("The provided DoFHandler does not have the same number of "
                            "components as the expected ones. For scalar fields there "
                            "should be 1 component."));

  // Set quadrature rule and FEValues to update the JxW values
  dealii::QGaussLobatto<dim> quadrature(degree + 1);
  dealii::FEValues<dim>      fe_values(dof_handler.get_fe(),
                                  quadrature,
                                  dealii::update_values | dealii::update_JxW_values);

  // Get the number of quadrature points
  const unsigned int num_quad_points = quadrature.size();

  // Create a value vector
  std::vector<number> quad_values(num_quad_points);

  // Loop over the cells provided by the DoFHandler
  number value = 0.0;
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Reinitialize the cell
          fe_values.reinit(cell);

          // Get the values
          fe_values.get_function_values(vector, quad_values);

          // Sum up the product of the JxW and values at each quadrature point to
          // compute the element integral.
          for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
            {
              value += quad_values[q_point] * fe_values.JxW(q_point);
            }
        }
    }

  integral_value = value;
}

template <int dim, int degree, typename number>
void
computeIntegral<dim, degree, number>::compute_integral(
  std::vector<number>           &integral_value,
  const dealii::DoFHandler<dim> &dof_handler,
  const VectorType              &vector) const
{
  [[maybe_unused]] const unsigned int expected_components = dim;
  Assert(dof_handler.get_fe().n_components() == expected_components,
         dealii::ExcMessage("The provided DoFHandler does not have the same number of "
                            "components as the expected ones. For vector fields there "
                            "should be dim components."));
  Assert(integral_value.size() == dim,
         dealii::ExcMessage("The provided `integral_value` must already be size dim"));

  // Set quadrature rule and FEValues to update the JxW values
  dealii::QGaussLobatto<dim> quadrature(degree + 1);
  dealii::FEValues<dim>      fe_values(dof_handler.get_fe(),
                                  quadrature,
                                  dealii::update_values | dealii::update_JxW_values);

  // Get the number of quadrature points
  const unsigned int num_quad_points = quadrature.size();

  // Create a value vector
  std::vector<dealii::Vector<number>> quad_values(num_quad_points,
                                                  dealii::Vector<number>(dim));

  // Loop over the cells and each lane in the vectorized array
  std::vector<number> value(dim, 0.0);
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
        {
          // Reinitialize the cell
          fe_values.reinit(cell);

          // Get the values
          fe_values.get_function_values(vector, quad_values);

          // Sum up the product of the JxW and values at each quadrature point to
          // compute the element integral.
          for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
            {
              for (unsigned int dimension = 0; dimension < dim; dimension++)
                {
                  value[dimension] +=
                    quad_values[q_point][dimension] * fe_values.JxW(q_point);
                }
            }
        }
    }

  for (unsigned int dimension = 0; dimension < dim; dimension++)
    {
      value[dimension] = dealii::Utilities::MPI::sum(value[dimension], MPI_COMM_WORLD);
    }

  integral_value = value;
}

INSTANTIATE_TRI_TEMPLATE(computeIntegral)

PRISMS_PF_END_NAMESPACE