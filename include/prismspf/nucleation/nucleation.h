// SPDX-FileCopyrightText: © 2025 PRISMS Center at the University of Michigan
// SPDX-License-Identifier: GNU Lesser General Public Version 2.1

#pragma once

#include <prismspf/core/types.h>

#include <prismspf/user_inputs/nucleation_parameters.h>

#include <prismspf/solvers/solver_context.h>

#include <prismspf/config.h>
#include <prismspf/nucleation/nucleus.h>

#include <algorithm>
#include <list>
#include <map>
#include <random>
#include <set>

PRISMS_PF_BEGIN_NAMESPACE

/**
 * @brief The class handles the stochastic nucleation
 * in PRISMS-PF.
 */
template <unsigned int dim, unsigned int degree, typename number>
class NucleationHandler
{
public:
  static void
  attempt_nucleation(const SolverContext<dim, degree, number> &solver_context,
                     std::vector<Nucleus<dim>>                &nuclei);

  static double
  calculate_number_of_events(const double &rate,
                             const double &delta_t,
                             const double &volume)
  {
    std::poisson_distribution<unsigned int> distribution(rate * delta_t * volume);
    return distribution(rng);
  }

private:
  // Add private member variables and functions as needed
  static std::uniform_real_distribution<double> uniform_unit_interval(-1.0, 1.0);
};

template <unsigned int dim, unsigned int degree, typename number>
inline static void
NucleationHandler<dim, degree, number>::attempt_nucleation(
  const SolverContext<dim, degree, number> &solver_context,
  std::vector<Nucleus<dim>>                &nuclei)
{
  NucleationParameters &nucleation_parameters =
    solver_context.user_inputs->get_nucleation_parameters();
  double delta_t =
    nucleation_parameters.get_nucleation_period() *
    solver_context.user_inputs->get_temporal_discretization().get_delta_t();
  auto &rng = solver_context.user_inputs->get_miscellaneous_parameters().rng;
  for (const auto &[index, variable] :
       solver_context.user_inputs->get_variable_attributes())
    {
      if (variable.is_nucleation_rate)
        {
          std::uniform_int_distribution<unsigned int> nucleating_index_dist(
            0,
            variable.nucleating_indices.size() - 1);
          std::list<Nucleus<dim>> new_nuclei;
          // Perform nucleation logic here
          // This is where you would check conditions and create nuclei
          for (const auto &cell : solver_context.get_triangulation_handler()
                                    .get_triangulation()
                                    .active_cell_iterators())
            {
              if (cell->is_locally_owned())
                {
                  // Grab the DoFHandler iterator
                  const auto dof_iterator = cell->as_dof_handler_iterator(
                    solver_context.get_dof_handler().get_dof_handler(index));

                  // Reinit the cell
                  fe_values.reinit(dof_iterator);
                  // Get the values for a scalar field
                  fe_values.get_function_values(
                    solver_context.get_solution_handler()
                      .get_solution_vector(index, DependencyType::Normal),
                    values);
                  for (unsigned int q_point = 0; q_point < num_quad_points; ++q_point)
                    {
                      nuc_rate +=
                        values[q_point] * fe_values.get_quadrature().weight(q_point);
                    }
                  unsigned int num_nuclei_in_cell =
                    nucleation_parameters.calculate_number_of_events(
                      nuc_rate,
                      delta_t,
                      solver_context.get_element_volume_container().get_volume(
                        dof_iterator));
                  for (unsigned int i = 0; i < num_nuclei_in_cell; ++i)
                    {
                      dealii::Point<dim> nucleus_location_unit_cell;
                      for (unsigned int d = 0; d < dim; ++d)
                        {
                          nucleus_location_unit_cell[d] = uniform_unit_interval(rng);
                        }
                      dealii::Point<dim> nucleus_location =
                        dealii::transform_unit_to_real_cell(cell,
                                                            nucleus_location_unit_cell);
                      double seed_time =
                        solver_context.user_inputs->get_temporal_discretization()
                          .get_time();
                      unsigned int seed_increment =
                        solver_context.user_inputs->get_temporal_discretization()
                          .get_increment();
                      unsigned int nucleating_index =
                        variable.nucleating_indices[nucleating_index_dist(rng)];

                      new_nuclei.emplace_back(nucleating_index,
                                              nucleus_location,
                                              seed_time,
                                              seed_increment);
                    }
                }
            }
          // Remove nuclei within their exclusion distance and add to nuclei list
          std::shuffle(new_nuclei.begin(), new_nuclei.end(), rng); // remove bias

          while (!new_nuclei.empty())
            {
              auto it = new_nuclei.begin();
              for (const auto &existing_nucleus : nuclei)
                {
                  if (it->location.distance(existing_nucleus.location) <
                        nucleation_parameters.get_exclusion_distance() ||
                      (it->field_index == existing_nucleus.field_index &&
                       it->location.distance(existing_nucleus.location) <
                         nucleation_parameters.get_same_field_exclusion_distance()))
                    {
                      new_nuclei.erase(it);
                      it = new_nuclei.end();
                      break;
                    }
                }
              if (it != new_nuclei.end())
                {
                  nuclei.push_back(*it);
                  new_nuclei.erase(it);
                }
            }
        }
    }
}

PRISMS_PF_END_NAMESPACE