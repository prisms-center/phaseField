// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <prismspf/core/field_container.h>

#include <prismspf/config.h>

#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief Compute the integral of a given field.
 */
template <unsigned int dim, unsigned int degree, typename number>
class Integrator
{
public:
  using VectorType = dealii::LinearAlgebra::distributed::Vector<number>;

  /**
   * @brief Constructor.
   */
  Integrator() = default;

  /**
   * @brief Compute the integral for a scalar field
   */
  void
  compute_integral(number                        &integral_value,
                   const dealii::DoFHandler<dim> &dof_handler,
                   const VectorType              &vector) const;

  /**
   * @brief Compute the integral for a vector field
   */
  void
  compute_integral(std::vector<number>           &integral_value,
                   const dealii::DoFHandler<dim> &dof_handler,
                   const VectorType              &vector) const;

  using TensorRank = FieldInfo::TensorRank;
  template <TensorRank Rank>
  using Value = std::conditional_t<Rank == TensorRank::Scalar,
                                   number,
                                   dealii::Tensor<int(Rank), dim, number>>;

  template <TensorRank rank>
  static Value<rank>
  integrate(const dealii::DoFHandler<dim> &dof_handler, const VectorType &vector)
  {
    [[maybe_unused]] constexpr unsigned int expected_components =
      dealii::Tensor<int(rank), dim>::n_independent_components;
    Assert(dof_handler.get_fe().n_components() == expected_components,
           dealii::ExcMessage("The provided DoFHandler does not have the same number of "
                              "components as the expected ones. For scalar fields there "
                              "should be 1 component."));

    // Update ghosts
    vector.update_ghost_values();

    // Set quadrature rule and FEValues to update the JxW values
    const dealii::QGaussLobatto<dim> quadrature(degree + 1);
    dealii::FEValues<dim>            fe_values(dof_handler.get_fe(),
                                    quadrature,
                                    dealii::update_values | dealii::update_JxW_values);

    // Get the number of quadrature points
    const unsigned int num_quad_points = quadrature.size();

    // Create a value vector
    std::vector<Value<rank>> quad_values(num_quad_points);

    // Loop over the cells provided by the DoFHandler
    Value<rank> value = 0;
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

    Value<rank> integral = dealii::Utilities::MPI::sum(value, MPI_COMM_WORLD);
    return integral;
  }
};

PRISMS_PF_END_NAMESPACE
