// SPDX-FileCopyrightText: © 2026 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <prismspf/core/system_wide.h>

#include <prismspf/utilities/assert.h>

#include <prismspf/config.h>

#include <functional>
#include <tuple>
#include <utility>
#include <vector>

PRISMS_PF_BEGIN_NAMESPACE

template <unsigned int dim, unsigned int degree, typename number>
class Integrator
{
public:
  using VectorType = dealii::LinearAlgebra::distributed::Vector<number>;

  /**
   * @brief Compute the integral of a given field.
   */
  template <int rank>
  static dealii::Vector<number>
  integrate(const dealii::DoFHandler<dim> &dof_handler, const auto &solution_vector)
  {
    constexpr unsigned int expected_components =
      dealii::Tensor<rank, dim>::n_independent_components;
    const unsigned int n_components = dof_handler.get_fe().n_components();

    DEBUG_ASSERT(n_components == expected_components,
                 "The provided DoFHandler does not have the same number of components as "
                 "the expected ones. For scalar fields there should be 1 component.");

    solution_vector.update_ghost_values();

    static dealii::FEValues<dim> fe_values(dof_handler.get_fe(),
                                           SystemWide<dim, degree>::quadrature,
                                           dealii::update_values |
                                             dealii::update_JxW_values);

    const unsigned int num_quad_points = SystemWide<dim, degree>::quadrature.size();

    std::vector<dealii::Vector<number>> quad_values(num_quad_points,
                                                    dealii::Vector<number>(
                                                      expected_components));

    dealii::Vector<number> value(expected_components);
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            fe_values.get_function_values(solution_vector, quad_values);

            // Sum up the product of the JxW and values at each quadrature point to
            // compute the element integral.
            for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
              {
                for (unsigned int component = 0; component < expected_components;
                     component++)
                  {
                    value[component] +=
                      quad_values[q_point][component] * fe_values.JxW(q_point);
                  }
              }
          }
      }

    for (unsigned int component = 0; component < expected_components; component++)
      {
        value[component] = dealii::Utilities::MPI::sum(value[component], MPI_COMM_WORLD);
      }
    return value;
  }

  /**
   * @brief Integrate function(u1(x),u2(x),...,un(x)) for a scalar fields.
   */
  template <typename Function, typename... SolutionVectors>
  static number
  integrate_function(const dealii::DoFHandler<dim> &dof_handler,
                     Function                     &&function,
                     const SolutionVectors &...solution_vectors)
  {
    // compile time check
    static_assert(sizeof...(SolutionVectors) > 0,
                  "At least one solution vector must be supplied.");

    DEBUG_ASSERT(dof_handler.get_fe().n_components() == 1,
                 "The input must be a scalar field.");

    (solution_vectors.update_ghost_values(), ...);

    const auto &quadrature = SystemWide<dim, degree>::quadrature;

    static dealii::FEValues<dim> fe_values(dof_handler.get_fe(),
                                           quadrature,
                                           dealii::update_values |
                                             dealii::update_JxW_values);

    const unsigned int num_quad_points = quadrature.size();

    const auto solutions = std::forward_as_tuple(solution_vectors...);

    auto quad_values = std::make_tuple(
      (static_cast<void>(solution_vectors), std::vector<number>(num_quad_points))...);

    number value = number {};

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);

            [&]<std::size_t... indices>(std::index_sequence<indices...>)
            {
              (fe_values.get_function_values(std::get<indices>(solutions),
                                             std::get<indices>(quad_values)),
               ...);
            }(std::index_sequence_for<SolutionVectors...> {});

            for (unsigned int q = 0; q < num_quad_points; ++q)
              {
                const number fvalue = std::apply(
                  [&](const auto &...values) -> number
                  {
                    return std::invoke(function, values[q]...);
                  },
                  quad_values);

                value += fvalue * fe_values.JxW(q);
              }
          }
      }

    return dealii::Utilities::MPI::sum(value, MPI_COMM_WORLD);
  }
};

PRISMS_PF_END_NAMESPACE
