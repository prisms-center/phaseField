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

    dealii::FEValues<dim> fe_values(dof_handler.get_fe(),
                                    SystemWide<dim, degree>::quadrature,
                                    dealii::update_values | dealii::update_JxW_values);

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
   * @brief Integrate function(u(x)) for a scalar field u.
   */
  template <typename Function, typename SolutionVector>
  static number
  integrate_function(const dealii::DoFHandler<dim> &dof_handler,
                     const SolutionVector          &solution_vector,
                     Function                     &&function)
  {
    DEBUG_ASSERT(dof_handler.get_fe().n_components() == 1,
                 "The input must be a scalar field.");

    solution_vector.update_ghost_values();

    dealii::FEValues<dim> fe_values(dof_handler.get_fe(),
                                    SystemWide<dim, degree>::quadrature,
                                    dealii::update_values | dealii::update_JxW_values);

    const unsigned int num_quad_points = SystemWide<dim, degree>::quadrature.size();

    std::vector<number> quad_values(num_quad_points);

    number value = number {};

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            fe_values.get_function_values(solution_vector, quad_values);

            for (unsigned int q = 0; q < num_quad_points; ++q)
              {
                value += std::invoke(function, quad_values[q]) * fe_values.JxW(q);
              }
          }
      }

    return dealii::Utilities::MPI::sum(value, MPI_COMM_WORLD);
  }

  /**
   * @brief Integrate function(u(x), v(x)) for two scalar fields.
   *
   * Both DoFHandlers must be defined on the same triangulation and have
   * corresponding active cells.
   */
  template <typename Function, typename SolutionVector1, typename SolutionVector2>
  static number
  integrate_function(const dealii::DoFHandler<dim> &dof_handler1,
                     const SolutionVector1         &solution1,
                     const dealii::DoFHandler<dim> &dof_handler2,
                     const SolutionVector2         &solution2,
                     Function                     &&function)
  {
    DEBUG_ASSERT(dof_handler1.get_fe().n_components() == 1,
                 "The first input must be a scalar field.");

    DEBUG_ASSERT(dof_handler2.get_fe().n_components() == 1,
                 "The second input must be a scalar field.");

    solution1.update_ghost_values();
    solution2.update_ghost_values();

    const auto &quadrature = SystemWide<dim, degree>::quadrature;

    dealii::FEValues<dim> fe_values1(dof_handler1.get_fe(),
                                     quadrature,
                                     dealii::update_values | dealii::update_JxW_values);

    dealii::FEValues<dim> fe_values2(dof_handler2.get_fe(),
                                     quadrature,
                                     dealii::update_values);

    const unsigned int num_quad_points = quadrature.size();

    std::vector<number> quad_values1(num_quad_points);
    std::vector<number> quad_values2(num_quad_points);

    number value = number {};

    auto cell1 = dof_handler1.begin_active();
    auto cell2 = dof_handler2.begin_active();

    for (; cell1 != dof_handler1.end() && cell2 != dof_handler2.end(); ++cell1, ++cell2)
      {
        DEBUG_ASSERT(cell1->id() == cell2->id(),
                     "The DoFHandlers do not have corresponding active cells.");

        if (cell1->is_locally_owned())
          {
            fe_values1.reinit(cell1);
            fe_values2.reinit(cell2);

            fe_values1.get_function_values(solution1, quad_values1);
            fe_values2.get_function_values(solution2, quad_values2);

            for (unsigned int q = 0; q < num_quad_points; ++q)
              {
                value += std::invoke(function, quad_values1[q], quad_values2[q]) *
                         fe_values1.JxW(q);
              }
          }
      }

    DEBUG_ASSERT(cell1 == dof_handler1.end() && cell2 == dof_handler2.end(),
                 "The DoFHandlers have different numbers of active cells.");

    return dealii::Utilities::MPI::sum(value, MPI_COMM_WORLD);
  }
};

PRISMS_PF_END_NAMESPACE
